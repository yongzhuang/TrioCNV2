package cn.edu.hit.triocnv.readdepth;

import cn.edu.hit.triocnv.readdepth.*;
import cn.edu.hit.triocnv.util.CNVRecord;
import cn.edu.hit.triocnv.util.Interval;
import cn.edu.hit.triocnv.util.IntervalReader;
import cn.edu.hit.triocnv.util.PEDReader;
import cn.edu.hit.triocnv.util.SVRecord;
import cn.edu.hit.triocnv.util.SVType;
import cn.edu.hit.triocnv.util.Trio;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;
import rcaller.RCaller;
import rcaller.RCode;

/**
 *
 * @author Yongzhuang Liu
 */
public class TrioCalling {

	private static org.apache.log4j.Logger logger = org.apache.log4j.Logger.getLogger(TrioCalling.class);
	private String inputFile;
	private String outputFile;
	private String pedFile;
	private double minMappability;
	int minDistance;
	int size = 500;
	private double outlier;
	double p;
	double e;
	double a = 0.0009;
	private String mappabilityFile;
	int gcBinSize;
	String[] samples;
	private String intervalFile;
	private static final int MIN_GC_BIN_SIZE = 1000;
	private static final int MIN_GC = 30;
	private static final int MAX_GC = 70;
	private static final int MAX_CONTIG_LENGTH = 10000000;

	public TrioCalling(String inputFile, String outputFile, String pedFile, String mappabilityFile,
			double minMappability, int minDistance, double transition_prob, double outlier, double e, int gcBinSize) {
		this.inputFile = inputFile;
		this.outputFile = outputFile;
		this.pedFile = pedFile;
		this.mappabilityFile = mappabilityFile;
		this.minMappability = minMappability;
		this.minDistance = minDistance;
		this.p = transition_prob;
		this.outlier = outlier;
		this.e = e;
		this.gcBinSize = gcBinSize;
	}

	private String getSexOfOffspring() {
		Trio trio = (new PEDReader(pedFile).getTrios()).get(0);
		if (trio.getOffspring().getSex() == 1) {
			return "M";
		} else if (trio.getOffspring().getSex() == 2) {
			return "F";
		} else {
			return "N";
		}
	}

	public void runMultiThreads(int numOfThreads) throws IOException {

		logger.info("Checking R environment and packages ......");
		try {
			Process pid = Runtime.getRuntime().exec("which Rscript");
			BufferedReader runTimeReader = new BufferedReader(new InputStreamReader(pid.getInputStream()));
			String R_HOME = runTimeReader.readLine().trim();
			if (R_HOME.equals("")) {
				logger.error("Rscript exectuable is not set in the PATH environment variable!");
				return;
			}
			if (!checkInstalledPackages(R_HOME)) {
				return;
			}
		} catch (Exception e) {
			e.printStackTrace();
		}

		logger.info("Loading data ......");
		List<Observation> observationList = getObservationList(inputFile);
		List<List<Observation>> observationListByChrom = getObservationListByChrom(observationList);
		if (this.intervalFile != null) {
			IntervalReader intervalReader = new IntervalReader(intervalFile);
			List<Interval> intervalList = intervalReader.getIntervalList();
			observationListByChrom = excludeRegions(observationListByChrom, intervalList);
		}

		List<Observation> observationListforParameterEstimation = new ArrayList();
		Pattern pattern = Pattern.compile("^(chr)?([1-9]|1[0-9]|2[0-2])$");
		for (int i = 0; i < observationListByChrom.size(); i++) {
			List<Observation> tmpObservationList = observationListByChrom.get(i);
			tmpObservationList = getFilteredObservationList(tmpObservationList, MIN_GC, MAX_GC, minMappability);
			observationListByChrom.set(i, tmpObservationList);
			String chrom = tmpObservationList.get(0).getChrom();
			Matcher matcher = pattern.matcher(chrom);
			if (matcher.matches()) {
				observationListforParameterEstimation.addAll(tmpObservationList);
			}
		}

		logger.info("Estimating parameters ......");
		NBModel[] nbModels = getNBModels(observationListByChrom.get(0));

		logger.info("Segmentation ......");
		double[] pi = getPI();
		double[][][][][][] transition = (new TransitionProbability(p, a, a)).getTransitionMatrix();
		File file = new File(outputFile);
		if (file.exists()) {
			try {
				file.delete();
				file.createNewFile();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
		ThreadPoolExecutor threadPool = (ThreadPoolExecutor) Executors.newFixedThreadPool(numOfThreads);
		List<List<Observation>> observationListByContig = getObservationListByContig(observationListByChrom,
				minDistance, MAX_CONTIG_LENGTH);
		for (int i = 0; i < observationListByContig.size(); i++) {
			List<Observation> observations = observationListByContig.get(i);
			String chrom = observations.get(0).getChrom();
			Matcher matcher = pattern.matcher(chrom);
			double[][][] inheritance = null;
			if (!matcher.matches()) {
				break;
			}
			if (matcher.matches()) {
				inheritance = (new InheritanceMatrix(e, a)).getAutoMatrix();
			} else if ((chrom.equals("X") || chrom.equals("chrX")) && getSexOfOffspring().equals("F")) {
				inheritance = (new InheritanceMatrix(e, a)).getFemaleXMatrix();
			} else if ((chrom.equals("X") || chrom.equals("chrX")) && getSexOfOffspring().equals("M")) {
				inheritance = (new InheritanceMatrix(e, a)).getMaleXMatrix();
			}
			Emission emission = new Emission(nbModels);
			TrioViterbi trioViterbi = (new TrioViterbi(observations, pi, transition, inheritance, emission));
			threadPool.execute(new SingleThreadCalling(trioViterbi, outputFile, minDistance));
		}
		threadPool.shutdown();
		try {
			threadPool.awaitTermination(30, TimeUnit.DAYS);
			boolean loop = true;
			long memory = 0;
			long maxTotalMemory = 0;
			do {
				long totalMemory = Runtime.getRuntime().totalMemory();
				long freeMemory = Runtime.getRuntime().freeMemory();
				long tmp = (totalMemory - freeMemory) / (1024 * 1024);
				if (tmp > memory) {
					memory = tmp;
				}
				if (totalMemory > maxTotalMemory) {
					maxTotalMemory = totalMemory;
				}
				loop = !threadPool.awaitTermination(1, TimeUnit.MINUTES);
			} while (loop);
			threadPool.awaitTermination(30, TimeUnit.DAYS);
		} catch (InterruptedException ex) {
			ex.printStackTrace();
		}
		printCNVRecord();
	}

	private void printCNVRecord() throws IOException {
		List<SVRecord> records = new ArrayList();
		try {
			BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(outputFile)));
			String line;
			while ((line = bufferedReader.readLine()) != null && line.trim().length() > 0) {
				String[] items = line.split("\t");
				String chrom = items[0];
				int start = Integer.parseInt(items[1]);
				int end = Integer.parseInt(items[2]);
				int[] states = new int[] { Integer.parseInt(items[3]), Integer.parseInt(items[4]),
						Integer.parseInt(items[5]) };
				double mappability = Double.parseDouble(items[6]);
				SVType[] typeArray=new SVType[3];
				typeArray[0]=parseSVType(states[0]);
				typeArray[1]=parseSVType(states[1]);
				typeArray[2]=parseSVType(states[2]);
				records.add(new SVRecord(chrom, start, end, typeArray, mappability));
			}
			bufferedReader.close();
		} catch (IOException ex) {
			ex.printStackTrace();
		}
		Collections.sort(records);
		String[] samples = getSampleName();
		try {
			PrintWriter writer = new PrintWriter(outputFile);
			writer.write("#CHROM" + "\t" + "START" + "\t" + "END" + "\t" + samples[0] + "\t" + samples[1] + "\t"
					+ samples[2] + "\t" + "Mappability" + "\t" +"Evidence"+ "\n");
			for (SVRecord record : records) {
				double mappability = getMappability(mappabilityFile, record.getChrom(), record.getStart(),
						record.getEnd());
				SVType[] types = record.getSVTypeArray();
				boolean refState = types[0] == SVType.REFERENCE && types[1] == SVType.REFERENCE  && types[2] == SVType.REFERENCE ;
				if (!refState) {
					writer.write(record.getChrom()+"\t"+record.getStart()+"\t"+record.getEnd()+"\t"+
				    types[0]+"\t"+types[1]+"\t"+types[2]+"\t"+mappability+"\t"+"RD"+"\n");
				}
			}
			writer.close();
		} catch (FileNotFoundException ex) {
			ex.printStackTrace();
		}

	}

	public double[] getPI() {
		double[] pi = new double[] { 0.0001, 0.0001, 0.00096, 0.0001, 0.0001 };
		return pi;
	}

	private String[] getSampleName() throws IOException {
		String[] samples = new String[3];
		BufferedReader rdFileReader = new BufferedReader(new FileReader(inputFile));
		String rdLine = rdFileReader.readLine();
		String[] record = rdLine.split("\t");
		samples[0] = record[3];
		samples[1] = record[4];
		samples[2] = record[5];
		rdFileReader.close();
		return samples;
	}

	public List<Observation> getObservationList(String readDepthFile) throws IOException {
		List<Observation> observationList = new ArrayList();
		BufferedReader rdFileReader = new BufferedReader(new FileReader(readDepthFile));
		String rdLine = rdFileReader.readLine();
		while ((rdLine = rdFileReader.readLine()) != null) {
			String[] record = rdLine.split("\t");
			String chrom = record[0];
			int start = Integer.parseInt(record[1]);
			int end = Integer.parseInt(record[2]);
			int[] rd = new int[] { Integer.parseInt(record[3]), Integer.parseInt(record[4]),
					Integer.parseInt(record[5]) };
			int gc = Integer.parseInt(record[6]);
			double mappability = Double.parseDouble(record[7]);
			Observation observation = new Observation(chrom, start, end, gc, mappability, rd);
			observationList.add(observation);
		}
		rdFileReader.close();
		return observationList;
	}

	private List<List<Observation>> getObservationListByContig(List<List<Observation>> observationListByChrom,
			int distance, int maxContigLength) {
		List<List<Observation>> observationListByContig = new ArrayList();
		for (List<Observation> observationList : observationListByChrom) {
			List<Observation> tmpObservationList = new ArrayList();
			tmpObservationList.add(observationList.get(0));
			for (int i = 1; i < observationList.size(); i++) {
				if (observationList.get(i).getEnd() - observationList.get(i - 1).getStart() + 1 > distance
						|| tmpObservationList.get(tmpObservationList.size() - 1).getEnd()
								- tmpObservationList.get(0).getStart() + 1 > maxContigLength) {
					observationListByContig.add(tmpObservationList);
					tmpObservationList = new ArrayList();
					tmpObservationList.add(observationList.get(i));
				} else {
					tmpObservationList.add(observationList.get(i));
				}
			}
			observationListByContig.add(tmpObservationList);
		}
		return observationListByContig;
	}

	private List<List<Observation>> getObservationListByChrom(List<Observation> observationList) {
		List<List<Observation>> observations = new ArrayList();
		List<Observation> tmp = null;
		String lastChrom = null;
		int lastEnd = 0;
		for (Observation observation : observationList) {
			String chrom = observation.getChrom();
			if (lastChrom == null) {
				tmp = new ArrayList();
				lastChrom = chrom;
			}
			if (lastChrom != null && !lastChrom.equals(chrom)) {
				observations.add(tmp);
				tmp = new ArrayList();
				lastChrom = chrom;
			}
			tmp.add(observation);
		}
		observations.add(tmp);
		return observations;
	}

	private List<Interval> getIntervalListOfChrom(List<Interval> intervalList, String chrom) {
		List<Interval> newIntervalList = new ArrayList();
		for (Interval interval : intervalList) {
			if (interval.getChrom().equals(chrom)) {
				newIntervalList.add(interval);
			}
		}
		return newIntervalList;
	}

	private List<List<Observation>> excludeRegions(List<List<Observation>> observationList,
			List<Interval> intervalList) {
		List<List<Observation>> newObservationList = new ArrayList();
		for (List<Observation> list : observationList) {
			List<Observation> currentList = new ArrayList();
			String chrom = list.get(0).getChrom();
			List<Interval> excludeList = getIntervalListOfChrom(intervalList, chrom);
			int i = 0;
			int j = 0;
			for (Interval interval : excludeList) {
				while (i < list.size() && list.get(i).getEnd() < interval.getStart()) {
					currentList.add(list.get(i));
					i++;
				}
				while (i < list.size() && list.get(i).getStart() >= interval.getStart()
						&& list.get(i).getStart() <= interval.getEnd()) {
					i++;
				}
			}
			while (i < list.size()) {
				currentList.add(list.get(i));
				i++;
			}
			newObservationList.add(currentList);
		}
		return newObservationList;
	}

	private List<Observation> getFilteredObservationList(List<Observation> observations, int minGC, int maxGC,
			double minMappability) {
		List<Observation> newObservations = new ArrayList();
		for (int j = 0; j < observations.size(); j++) {
			Observation tmp = observations.get(j);
			if (tmp.getGC() >= minGC && tmp.getGC() <= maxGC && tmp.getMappability() >= minMappability) {
				newObservations.add(tmp);
			}
		}
		return newObservations;
	}

	private NBModel[] getNBModels(List<Observation> observationList) {
		int size = observationList.size();
		int[] fatherRD = new int[size];
		int[] motherRD = new int[size];
		int[] offspringRD = new int[size];
		int[] gc = new int[size];
		double[] mappability = new double[size];
		for (int i = 0; i < observationList.size(); i++) {
			Observation observation = observationList.get(i);
			fatherRD[i] = observation.getRD()[0];
			motherRD[i] = observation.getRD()[1];
			offspringRD[i] = observation.getRD()[2];
			gc[i] = observation.getGC();
			mappability[i] = observation.getMappability();
		}
		try {
			NBModel[] nbModels = new NBModel[3];
			Process pid = Runtime.getRuntime().exec("which Rscript");
			BufferedReader runTimeReader = new BufferedReader(new InputStreamReader(pid.getInputStream()));
			String R_HOME = runTimeReader.readLine().trim();
			runTimeReader.close();
			double[] outliers = new double[] { outlier / 2, 1 - outlier / 2 };
			int[] binsize = new int[] { MIN_GC_BIN_SIZE };
			for (int k = 0; k < 3; k++) {
				int[] sample = new int[] { k + 1 };
				RCaller caller = new RCaller();
				RCode code = new RCode();
				caller.setRscriptExecutable(R_HOME);
				caller.setRCode(code);
				code.addIntArray("sample", sample);
				code.addIntArray("fatherRD", fatherRD);
				code.addIntArray("motherRD", motherRD);
				code.addIntArray("offspringRD", motherRD);
				code.addIntArray("gc", gc);
				code.addDoubleArray("mappability", mappability);
				code.addDoubleArray("outliers", outliers);
				code.addIntArray("size", binsize);
				code.addRCode("library(MASS)");
				code.addRCode("data<-data.frame(fatherRD,motherRD,offspringRD,gc,mappability)");
				code.addRCode("data<-split(data[,c(sample[1],5)],data$gc)");
				code.addRCode("theta<-rep(0,101)");
				code.addRCode("coef1<-rep(0,101)");
				code.addRCode("coef2<-rep(0,101)");
				code.addRCode("remove_outliers <- function(data, na.rm = TRUE, ...) {");
				code.addRCode("range <- quantile(data[,1], prob=c(outliers[1], outliers[2]),na.rm = TRUE)");
				code.addRCode("new <- data");
				code.addRCode("new[which(new[,1] >= range[1] &  new[,1] <= range[2]),]");
				code.addRCode("}");
				code.addRCode("for(i in 1:length(data)){");
				code.addRCode("bin<-as.data.frame(unname(data[i]))");
				code.addRCode("if(nrow(bin)>size[1]){");
				code.addRCode("bin<-remove_outliers(bin)");
				code.addRCode("index<-sample(1:nrow(bin), size = min(nrow(bin),size[1]))");
				code.addRCode("count<-bin[index,1]");
				code.addRCode("mappability<-bin[index,2]");
				code.addRCode("fit<-glm.nb(count~mappability)");
				code.addRCode("theta[as.numeric(names(data[i]))+1]<-fit$theta");
				code.addRCode("coef1[as.numeric(names(data[i]))+1]<-unname(fit$coef[1])");
				code.addRCode("coef2[as.numeric(names(data[i]))+1]<-unname(fit$coef[2])");
				code.addRCode("}");
				code.addRCode("}");
				code.addRCode("result<-list(theta=theta,coef1=coef1, coef2=coef2)");
				caller.runAndReturnResult("result");
				double[] theta = caller.getParser().getAsDoubleArray("theta");
				double[] coef1 = caller.getParser().getAsDoubleArray("coef1");
				double[] coef2 = caller.getParser().getAsDoubleArray("coef2");
				nbModels[k] = new NBModel(theta, coef1, coef2);
			}
			return nbModels;
		} catch (IOException ex) {
			ex.printStackTrace();
			return null;
		}
	}

	private String toCNVType(int state) {
		String type = null;
		switch (state) {
		case 0:
			type = "DEL";
			break;
		case 1:
			type = "DEL";
			break;
		case 2:
			type = "REF";
			break;
		case 3:
			type = "DUP";
			break;
		case 4:
			type = "DUP";
			break;
		}
		return type;
	}

	private double getMappability(String mappabilityFile, String chrom, int start, int end) {
		if (!chrom.startsWith("chr")) {
			chrom = "chr" + chrom;
		}
		try {
			BBFileReader reader = new BBFileReader(mappabilityFile);
			BigWigIterator iterator = reader.getBigWigIterator(chrom, start - 1, chrom, end, false);
			double total = 0;
			while (iterator.hasNext()) {
				WigItem item = iterator.next();
				double value = item.getWigValue();
				int startBase = item.getStartBase();
				int endBase = item.getEndBase();
				if (startBase < start - 1 && endBase > end) {
					total += value * (end - start + 1);
				} else if (startBase < start - 1 && endBase <= end) {
					total += value * (endBase - start + 1);
				} else if (startBase >= start - 1 && endBase > end) {
					total += value * (end - startBase);
				} else {
					total += value * (endBase - startBase);
				}
			}
			reader.close();
			return total / (end - start + 1);
		} catch (IOException ex) {
			ex.printStackTrace();
			return 0;
		}
	}

	private boolean checkInstalledPackages(String R_HOME) {
		boolean tag = true;
		RCaller caller = new RCaller();
		caller.setRscriptExecutable(R_HOME);
		caller.cleanRCode();
		RCode code = new RCode();
		String[] packages = new String[] { "Runiversal", "MASS" };
		code.addStringArray("packages", packages);
		code.addRCode("label<-c(1, 1)");
		code.addRCode("for(i in 1:2){");
		code.addRCode("if(!require(packages[i], character.only=TRUE)){");
		code.addRCode("label[i]=0");
		code.addRCode("}");
		code.addRCode("}");
		caller.setRCode(code);
		caller.runAndReturnResult("label");
		int[] label = caller.getParser().getAsIntArray("label");
		for (int i = 0; i < packages.length; i++) {
			if (label[i] == 0) {
				logger.error(packages[i] + " is not installed in R environment!");
				tag = false;
			}
		}
		return tag;
	}
	
	private SVType parseSVType(int state) {
		if (state<2) {
			return SVType.DELETION;
		}
		else if (state==2) {
			return SVType.REFERENCE;
		}
		else{
			return SVType.DUPLICATION;
		}
	}
}
