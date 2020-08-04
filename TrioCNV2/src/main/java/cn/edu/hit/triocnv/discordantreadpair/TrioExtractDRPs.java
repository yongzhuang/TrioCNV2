package cn.edu.hit.triocnv.discordantreadpair;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;

import cn.edu.hit.triocnv.readdepth.Preprocessing;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;

/**
 *
 * @author Yongzhuang Liu
 */
public class TrioExtractDRPs {

	private static Logger logger = Logger.getLogger(TrioExtractDRPs.class);
	private String referenceFile;
	private String outputFolder;
	private List<String> bamFileList;
	private static final int MAX_DISTANCE = 10000000;

	public TrioExtractDRPs(String referenceFile, List<String> bamFileList, String outputFolder) {
		this.referenceFile = referenceFile;
		this.outputFolder = outputFolder;
		this.bamFileList = bamFileList;
	}

	public void extract(Map<String, int[]> parameterMap, int deviation) throws IOException {
		ReferenceSource referenceSource = new ReferenceSource(new File(referenceFile));
		SamReader samReader = SamReaderFactory.makeDefault().referenceSource(referenceSource)
				.open(new File(bamFileList.get(0)));
		SAMSequenceDictionary samSequenceDictionary = samReader.getFileHeader().getSequenceDictionary();
		Pattern pattern = Pattern.compile("^(chr)?([1-9]|1[0-9]|2[0-2])$");
		for (int i = 0; i < samSequenceDictionary.getSequences().size(); i++) {
			SAMSequenceRecord samSequenceRecord = samSequenceDictionary.getSequence(i);
			String sequenceName = samSequenceRecord.getSequenceName();
			Matcher matcher = pattern.matcher(sequenceName);
			if (matcher.matches()) {
				logger.info("Chromosome " + sequenceName + " is processing ... ...");
				List<ReadPair> delReadPairList = new ArrayList();
				List<ReadPair> dupReadPairList = new ArrayList();
				for (String bamFile : bamFileList) {
					samReader = SamReaderFactory.makeDefault().referenceSource(referenceSource)
							.open(new File(bamFileList.get(0)));
					String sample = samReader.getFileHeader().getReadGroups().get(0).getSample();
					int[] parameters = parameterMap.get(sample);
					List<List<ReadPair>> drpList = extractDiscordantReadPairs(bamFile, sequenceName, 1,
							samSequenceRecord.getSequenceLength(), parameters[0], parameters[1], deviation);
					delReadPairList.addAll(drpList.get(0));
					dupReadPairList.addAll(drpList.get(1));
				}
				Collections.sort(delReadPairList);
				Collections.sort(dupReadPairList);
				File outputFolder = new File(this.outputFolder);
				if (outputFolder.exists()) {
					if (!outputFolder.isDirectory()) {
						return;
					}
				} else {
					outputFolder.mkdir();
				}
				FileWriter deletionWriter = new FileWriter(new File(outputFolder + "/" + sequenceName + ".DEL"));
				for (ReadPair readPair : delReadPairList) {
					deletionWriter.write(readPair + "\n");
				}
				deletionWriter.close();
				FileWriter duplicationWriter = new FileWriter(new File(outputFolder + "/" + sequenceName + ".DUP"));
				for (ReadPair readPair : dupReadPairList) {
					duplicationWriter.write(readPair + "\n");
				}
				duplicationWriter.close();
			}
		}
	}

	public List<List<ReadPair>> extractDiscordantReadPairs(String bamFile, String chrom, int start, int end, int median,
			int mad, int deviation) throws IOException {
		List<List<ReadPair>> drpList = new ArrayList();
		ReferenceSource referenceSource = new ReferenceSource(new File(referenceFile));
		SamReader samReader = SamReaderFactory.makeDefault().referenceSource(referenceSource).open(new File(bamFile));
		SAMRecordIterator samRecordIterator = samReader.query(chrom, start, end, false);
		HashMap<String, SAMRecord> delSAMRecordTable = new HashMap<String, SAMRecord>();
		HashMap<String, SAMRecord> dupSAMRecordTable = new HashMap<String, SAMRecord>();
		List<ReadPair> delReadPairList = new ArrayList();
		List<ReadPair> dupReadPairList = new ArrayList();
		while (samRecordIterator.hasNext()) {
			SAMRecord samRecord = samRecordIterator.next();
			if (!isCandidateDRP(samRecord)) {
				continue;
			}
			String readName = samRecord.getReadName();
			int insertSize = Math.abs(samRecord.getInferredInsertSize());
			SamPairUtil.PairOrientation orientation = SamPairUtil.getPairOrientation(samRecord);
			if (orientation == SamPairUtil.PairOrientation.FR) {
				if (insertSize > median + mad * deviation) {
					if (delSAMRecordTable.containsKey(readName)) {
						SAMRecord leftSAMRecord = delSAMRecordTable.get(readName);
						String sample = samRecord.getReadGroup().getSample();
						int leftStart = leftSAMRecord.getAlignmentStart();
						int leftEnd = leftSAMRecord.getAlignmentEnd();
						int rightStart = samRecord.getAlignmentStart();
						int rightEnd = samRecord.getAlignmentEnd();
						delReadPairList.add(new ReadPair(sample, chrom, leftStart, leftEnd, rightStart, rightEnd));
						delSAMRecordTable.remove(readName);
					} else {
						delSAMRecordTable.put(readName, samRecord);
					}
				}
			}
			if (orientation == SamPairUtil.PairOrientation.RF) {
				if (dupSAMRecordTable.containsKey(readName)) {
					SAMRecord leftSAMRecord = dupSAMRecordTable.get(readName);
					String sample = samRecord.getReadGroup().getSample();
					int leftStart = leftSAMRecord.getAlignmentStart();
					int leftEnd = leftSAMRecord.getAlignmentEnd();
					int rightStart = samRecord.getAlignmentStart();
					int rightEnd = samRecord.getAlignmentEnd();
					dupReadPairList.add(new ReadPair(sample, chrom, leftStart, leftEnd, rightStart, rightEnd));
					dupSAMRecordTable.remove(readName);
				} else {
					dupSAMRecordTable.put(readName, samRecord);
				}
			}
		}
		samRecordIterator.close();
		samReader.close();
		Collections.sort(delReadPairList);
		Collections.sort(dupReadPairList);
		drpList.add(delReadPairList);
		drpList.add(dupReadPairList);
		return drpList;
	}

	private boolean isCandidateDRP(SAMRecord samRecord) {
		if (!samRecord.getReadPairedFlag() || samRecord.getReadUnmappedFlag() || samRecord.getDuplicateReadFlag()
				|| samRecord.getMateUnmappedFlag() || samRecord.isSecondaryOrSupplementary()
				|| samRecord.getReferenceIndex() != samRecord.getMateReferenceIndex()
				|| Math.abs(samRecord.getMateAlignmentStart() - samRecord.getAlignmentStart()) > MAX_DISTANCE
				|| samRecord.getMappingQuality() == 0) {
			return false;
		} else {
			return true;
		}
	}
}
