/**
 * 
 */
package cn.edu.hit.triocnv.breakpoint;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;

import cn.edu.hit.triocnv.discordantreadpair.EstimateInsertSizes;
import cn.edu.hit.triocnv.readdepth.Emission;
import cn.edu.hit.triocnv.readdepth.InheritanceMatrix;
import cn.edu.hit.triocnv.readdepth.Observation;
import cn.edu.hit.triocnv.readdepth.SingleThreadCalling;
import cn.edu.hit.triocnv.readdepth.TrioViterbi;
import cn.edu.hit.triocnv.util.PEDReader;
import cn.edu.hit.triocnv.util.SVRecord;
import cn.edu.hit.triocnv.util.SVType;
import cn.edu.hit.triocnv.util.Trio;

/**
 * @author Yongzhuang Liu
 *
 */
public class MultiThreadBreakPointCalling {

	private String referenceFile;
	private String bamList;
	private String pedFile;
	private String rdSVFile;
	private String drpSVFile;
	private String workDir;
	private int numOfThreads;
	private int size;
	private int deviation;

	public MultiThreadBreakPointCalling(String referenceFile, String bamList, String pedFile, String rdSVFile,
			String drpSVFile, String workDir, int size, int deviation, int numOfThreads) {
		this.referenceFile = referenceFile;
		this.bamList = bamList;
		this.pedFile = pedFile;
		this.rdSVFile = rdSVFile;
		this.drpSVFile = drpSVFile;
		this.workDir = workDir;
		this.size = size;
		this.deviation = deviation;
		this.numOfThreads = numOfThreads;
	}

	public void run() throws IOException, InterruptedException {
		List<SVRecord> drpSVRecordList = getSVRecordList(drpSVFile);
		List<SVRecord> rdSVRecordList = getSVRecordList(rdSVFile);
		List<SVRecord> mergedSVRecordList = merge(drpSVRecordList, rdSVRecordList, 0.5);
		List<String> bamFileList = new ArrayList();
		BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(bamList)));
		String line = null;
		while ((line = bufferedReader.readLine()) != null && line.trim().length() > 0) {
			bamFileList.add(line);
		}
		bufferedReader.close();
		String parameterFile = workDir + "/parameters.txt";
		EstimateInsertSizes estimateInsertSizes = new EstimateInsertSizes(referenceFile, bamFileList, parameterFile);
		estimateInsertSizes.estimate();

		ThreadPoolExecutor threadPool = (ThreadPoolExecutor) Executors.newFixedThreadPool(numOfThreads);
		List<List<SVRecord>> list = subListBySegment(mergedSVRecordList, numOfThreads);
		for (int i = 0; i < list.size(); i++) {
			threadPool.execute(new SingleThreadRefinement(list.get(i), referenceFile, bamFileList, pedFile,
					parameterFile, workDir, size, deviation));
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
		outputSVRecord();
	}

	private void outputSVRecord() throws IOException {
		List<SVRecord> resultSVRecordList = new ArrayList();
		BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(workDir + "/Temp_CNV.txt")));
		String line = null;
		while ((line = bufferedReader.readLine()) != null && line.trim().length() > 0) {
			if (!line.startsWith("#")) {
				String[] record = line.split("\t");
				String chrom = record[0];
				int start = Integer.parseInt(record[1]);
				int end = Integer.parseInt(record[2]);
				SVType[] svTypeArray = new SVType[3];
				svTypeArray[0] = SVType.valueOf(record[3]);
				svTypeArray[1] = SVType.valueOf(record[4]);
				svTypeArray[2] = SVType.valueOf(record[5]);
				String evidence = record[6];
				resultSVRecordList.add(new SVRecord(chrom, start, end, svTypeArray, evidence));
			}
		}
		bufferedReader.close();
		Collections.sort(resultSVRecordList);
		Trio trio = (new PEDReader(pedFile).getTrios()).get(0);
		String fatherID = trio.getFather().getIndividualID();
		String motherID = trio.getMother().getIndividualID();
		String offspringID = trio.getOffspring().getIndividualID();
		PrintWriter writer = new PrintWriter(workDir + "/CNV.txt");
		writer.write("#CHROM\tSTART\tEND\t" + fatherID + "\t" + motherID + "\t" + offspringID + "\tEvidence\n");
		for (SVRecord tmpSVRecord : resultSVRecordList) {
			writer.write(tmpSVRecord.getChrom() + "\t" + tmpSVRecord.getStart() + "\t" + tmpSVRecord.getEnd() + "\t"
					+ tmpSVRecord.getSVTypeArray()[0] + "\t" + tmpSVRecord.getSVTypeArray()[1] + "\t"
					+ tmpSVRecord.getSVTypeArray()[2] + "\t" + tmpSVRecord.getEvidence() + "\n");
		}
		writer.close();
	}

	private List<SVRecord> merge(List<SVRecord> svRecordList1, List<SVRecord> svRecordList2, double overlap) {
		String[] chroms = { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
				"18", "19", "20", "21", "22", "X", "Y" };
		List<SVRecord> svRecordList = new ArrayList();
		for (String chrom : chroms) {
			List<SVRecord> tmpSVRecordList1 = getSVRecordOfChrom(svRecordList1, chrom);
			List<SVRecord> tmpSVRecordList2 = getSVRecordOfChrom(svRecordList2, chrom);
			int i = 0;
			int j = 0;
			while (i < tmpSVRecordList1.size() && j < tmpSVRecordList2.size()) {
				if (tmpSVRecordList1.get(i).getStart() > tmpSVRecordList2.get(j).getEnd()) {
					svRecordList.add(tmpSVRecordList2.get(j++));
				} else if (tmpSVRecordList2.get(j).getStart() > tmpSVRecordList1.get(i).getEnd()) {
					svRecordList.add(tmpSVRecordList1.get(i++));
				} else {
					if (isOverlap(tmpSVRecordList1.get(i), tmpSVRecordList2.get(j), overlap)) {
						SVRecord svRecord = tmpSVRecordList1.get(i);
						svRecord.setEvidence("RD:RP");
						svRecordList.add(svRecord);
					} else {
						if (tmpSVRecordList1.get(i).getStart() >= tmpSVRecordList2.get(j).getStart()) {
							svRecordList.add(tmpSVRecordList2.get(j));
							svRecordList.add(tmpSVRecordList1.get(i));
						} else {
							svRecordList.add(tmpSVRecordList1.get(i));
							svRecordList.add(tmpSVRecordList2.get(j));
						}
					}
					i++;
					j++;
				}
			}
			while (i < tmpSVRecordList1.size()) {
				svRecordList.add(tmpSVRecordList1.get(i++));
			}
			while (j < tmpSVRecordList2.size()) {
				svRecordList.add(tmpSVRecordList2.get(j++));
			}
		}
		return svRecordList;
	}

	private boolean isOverlap(SVRecord svRecord1, SVRecord svRecord2, double percent) {
		if (!svRecord1.getChrom().equals(svRecord2.getChrom())) {
			return false;
		} else if (svRecord1.getStart() > svRecord2.getEnd() || svRecord2.getStart() > svRecord1.getEnd()) {
			return false;
		} else {
			double overlap = 0;
			if (svRecord1.getStart() >= svRecord2.getStart()) {
				if (svRecord1.getEnd() <= svRecord2.getEnd()) {
					overlap = svRecord1.getLength();
				} else {
					overlap = svRecord2.getEnd() - svRecord1.getStart() + 1;
				}
			} else {
				if (svRecord2.getEnd() <= svRecord1.getEnd()) {
					overlap = svRecord2.getLength();
				} else {
					overlap = svRecord1.getEnd() - svRecord2.getStart() + 1;
				}
			}
			if (overlap / (double) svRecord1.getLength() > percent
					&& overlap / (double) svRecord2.getLength() > percent) {
				return true;
			} else {
				return false;
			}
		}
	}

	private List<SVRecord> getSVRecordOfChrom(List<SVRecord> svRecordList, String chrom) {
		List<SVRecord> resultSVRecordList = new ArrayList();
		for (SVRecord svRecord : svRecordList) {
			if (svRecord.getChrom().equals(chrom) || svRecord.getChrom().equals("chr" + chrom)) {
				resultSVRecordList.add(svRecord);
			}
		}
		return resultSVRecordList;
	}

	public List<SVRecord> getSVRecordList(String svFile) throws IOException {
		List<SVRecord> svRecordList = new ArrayList();
		BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(svFile)));
		String line = null;
		while ((line = bufferedReader.readLine()) != null && line.trim().length() > 0) {
			if (!line.startsWith("#")) {
				String[] record = line.split("\t");
				String chrom = record[0];
				int start = Integer.parseInt(record[1]);
				int end = Integer.parseInt(record[2]);
				SVType[] svTypeArray = new SVType[3];
				svTypeArray[0] = SVType.valueOf(record[3]);
				svTypeArray[1] = SVType.valueOf(record[4]);
				svTypeArray[2] = SVType.valueOf(record[5]);
				String evidence = record[7];
				svRecordList.add(new SVRecord(chrom, start, end, svTypeArray, evidence));
			}
		}
		bufferedReader.close();
		return svRecordList;
	}

	private List<List<SVRecord>> subListBySegment(List<SVRecord> data, int segments) {
		List<List<SVRecord>> result = new ArrayList<>();
		int size = data.size();
		if (size > 0 && segments > 0) {
			int count = size / segments;
			List<SVRecord> cutList = null;
			for (int i = 0; i < segments; i++) {
				if (i == segments - 1) {
					cutList = data.subList(count * i, size);
				} else {
					cutList = data.subList(count * i, count * (i + 1));
				}
				result.add(cutList);
			}
		} else {
			result.add(data);
		}
		return result;
	}
}
