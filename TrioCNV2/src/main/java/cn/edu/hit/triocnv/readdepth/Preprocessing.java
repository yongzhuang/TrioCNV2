package cn.edu.hit.triocnv.readdepth;

import cn.edu.hit.triocnv.util.MultiSamRecordIterator;
import cn.edu.hit.triocnv.util.PEDReader;
import cn.edu.hit.triocnv.util.Trio;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.URI;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.log4j.Logger;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

/**
 *
 * @author Yongzhuang Liu
 */
public class Preprocessing {

	private static Logger logger = Logger.getLogger(Preprocessing.class);
	private String referenceSequenceFile;
	private String bamsFile;
	private String pedFile;
	private String mappabilityFile;
	private String outputFile;
	private int windowSize;
	private int minMappingQuality;

	public Preprocessing(String referenceSequenceFile, String bamsFile, String pedFile, String mappabilityFile,
			String outputFile, int windowSize, int minMappingQuality) {
		this.referenceSequenceFile = referenceSequenceFile;
		this.bamsFile = bamsFile;
		this.pedFile = pedFile;
		this.mappabilityFile = mappabilityFile;
		this.outputFile = outputFile;
		this.windowSize = windowSize;
		this.minMappingQuality = minMappingQuality;
	}

	public void process() {
		try {
			Trio trio = (new PEDReader(pedFile).getTrios()).get(0);
			String[] samples = new String[3];
			samples[0] = trio.getFather().getIndividualID();
			samples[1] = trio.getMother().getIndividualID();
			samples[2] = trio.getOffspring().getIndividualID();
			Map<String, Integer> index = new HashMap();
			for (int i = 0; i < 3; i++) {
				index.put(samples[i], i);
			}
			ReferenceSource referenceSource = new ReferenceSource(new File(referenceSequenceFile));
			List<SamReader> readers = new ArrayList();
			BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(bamsFile)));
			String line = null;
			while ((line = bufferedReader.readLine()) != null && line.trim().length() > 0) {
				SamReader samReader = SamReaderFactory.makeDefault().referenceSource(referenceSource)
						.open(new File(line));
				readers.add(samReader);
			}
			bufferedReader.close();
			PrintWriter writer = new PrintWriter(new File(outputFile));
			writer.write("CHROM\tSTART\tEND");
			for (int i = 0; i < samples.length; i++) {
				writer.write("\t" + samples[i]);
			}
			writer.write("\tGC\tMappability\n");
			ReferenceSequenceFileWalker walker = new ReferenceSequenceFileWalker(new File(referenceSequenceFile));
			int[] RD = new int[readers.size()];
			MultiSamRecordIterator iterators = new MultiSamRecordIterator(readers);
			int bin = -1;
			String chrom = null;
			byte[] bases = null;
			double[] mappabilities = null;
			while (iterators.hasNext()) {
				SAMRecord nextRead = iterators.next();
				if (nextRead.isSecondaryOrSupplementary() || nextRead.getReadUnmappedFlag()
						|| nextRead.getDuplicateReadFlag() || nextRead.getMappingQuality() < minMappingQuality) {
					continue;
				}
				String id = nextRead.getReadGroup().getSample();
				String referenceName = nextRead.getReferenceName();
				int referenceIndex = nextRead.getReferenceIndex();
				Pattern pattern = Pattern.compile("^(chr)?([1-9]|1[0-9]|2[0-2]|[X|Y])$");
				Matcher matcher = pattern.matcher(referenceName);
				if (!matcher.matches()) {
					break;
				}
				int start = nextRead.getAlignmentStart();
				int currentBIN = (start - 1) / windowSize;
				if (chrom == null || !referenceName.equals(chrom)) {
					bin = currentBIN;
					RD = new int[readers.size()];
					chrom = referenceName;
					bases = walker.get(referenceIndex).getBases();
					int length = bases.length;
					mappabilities = getMappability(mappabilityFile, chrom, length, windowSize);
					logger.info("Chromosome " + chrom + " is processing ... ...");
				}
				int GC;
				while (currentBIN != bin) {
					GC = (new GCCalculation(bases, bin * windowSize + 1, (bin + 1) * windowSize + 1)).getGCContent();
					double mappability = mappabilities[bin];
					writer.write(buildRecord(chrom, bin, windowSize, RD, GC, mappability) + "\n");
					bin = bin + 1;
					RD = new int[readers.size()];
				}
				RD[index.get(id)]++;
			}
			writer.close();
		} catch (IOException ex) {
			ex.printStackTrace();
		}
	}

	private String buildRecord(String chrom, int bin, int window, int[] RD, int GC, double mappability) {
		String result = chrom + "\t" + (bin * window + 1) + "\t" + ((bin + 1) * window);
		for (int i = 0; i < RD.length; i++) {
			result += "\t" + RD[i];
		}
		result += "\t" + GC + "\t" + mappability;
		return result;
	}

	private double[] getMappability(String mappabilityFile, String chrom, int length, int windowSize) {
		if (!chrom.startsWith("chr")) {
			chrom = "chr" + chrom;
		}
		try {
			List<Double> mappability = new ArrayList();
			BBFileReader reader = new BBFileReader(mappabilityFile);
			BigWigIterator iterator = reader.getBigWigIterator(chrom, 1, chrom, length, true);
			double[] result = new double[length / windowSize];
			int bin = 0;
			double tmp = 0;
			while (iterator.hasNext()) {
				WigItem item = iterator.next();
				int start = item.getStartBase();
				int end = item.getEndBase();
				double value = item.getWigValue();
				for (int i = start + 1; i <= end; i++) {
					int currentBin = (i - 1) / windowSize;
					if (currentBin != bin) {
						if (tmp > 0) {
							result[bin] = tmp / windowSize;
							tmp = 0;
						}
						bin = currentBin;
					}
					tmp += value;
				}
			}
			return result;
		} catch (IOException ex) {
			ex.printStackTrace();
			return null;
		}
	}
}
