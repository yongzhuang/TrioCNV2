package cn.edu.hit.triocnv.discordantreadpair;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Median;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;

/**
*
* @author Yongzhuang Liu
*/

public class EstimateInsertSizes {

	private String referenceFile;
	private List<String> bamFileList;
	private static final int NUM_READ_PAIRS = 1000000;
	private static final double OUTLIER = 0.05;
	private String outputFile;

	public EstimateInsertSizes(String referenceFile, List<String> bamFileList, String outputFile) {
		this.referenceFile = referenceFile;
		this.bamFileList = bamFileList;
		this.outputFile = outputFile;
	}

	public void estimate() throws IOException {
		FileWriter parameterWriter = new FileWriter(new File(outputFile));
		for (String bamFile : bamFileList) {
			ReferenceSource referenceSource = new ReferenceSource(new File(referenceFile));
			SamReader samReader = SamReaderFactory.makeDefault().referenceSource(referenceSource)
					.open(new File(bamFile));
			String sample = samReader.getFileHeader().getReadGroups().get(0).getSample();
			SAMRecordIterator samRecordIterator = samReader.iterator();
			double[] insertSizeArray = new double[NUM_READ_PAIRS];
			int i = 0;
			while (samRecordIterator.hasNext()) {
				SAMRecord samRecord = samRecordIterator.next();
				if (!samRecord.getReadPairedFlag() || !samRecord.getFirstOfPairFlag()
						|| samRecord.getDuplicateReadFlag() || samRecord.getReadUnmappedFlag()
						|| samRecord.getMateUnmappedFlag() || samRecord.getSupplementaryAlignmentFlag()
						|| samRecord.getAttribute("XA") != null) {
					continue;
				}
				SamPairUtil.PairOrientation orientation = SamPairUtil.getPairOrientation(samRecord);
				if (orientation == SamPairUtil.PairOrientation.FR) {
					int insertSize = Math.abs(samRecord.getInferredInsertSize());
					if (i < NUM_READ_PAIRS) {
						insertSizeArray[i] = insertSize;
						i++;
					} else {
						break;
					}
				}
			}
			samRecordIterator.close();
			samReader.close();
			double[] newInsertSizeArray = getOutliersRemovedArray(insertSizeArray, OUTLIER, 1 - OUTLIER);
			double median = (new Median()).evaluate(newInsertSizeArray);
			double standardDeviation = (new StandardDeviation()).evaluate(newInsertSizeArray);
			parameterWriter.write(sample + "\t" + median + "\t" + standardDeviation + "\n");
		}
		parameterWriter.close();
	}

	private double[] getOutliersRemovedArray(double[] data, double minP, double maxP) {
		int n = data.length;
		Arrays.sort(data);
		int minIndex = (int) (minP * n) - 1;
		int maxIndex = (int) (maxP * n) - 1;
		double[] data2 = new double[maxIndex - minIndex + 1];
		for (int i = minIndex; i <= maxIndex; i++) {
			data2[i - minIndex] = data[i];
		}
		return data2;
	}
}
