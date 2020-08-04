package cn.edu.hit.triocnv.breakpoint;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import cn.edu.hit.triocnv.util.SVRecord;
import cn.edu.hit.triocnv.util.SVType;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqEncoder;
import htsjdk.samtools.fastq.FastqRecord;

/**
 *
 * @author Yongzhuang Liu
 */

public class SoftClipCalling {

	private String referenceFile;
	private List<String> bamFileList;

	public SoftClipCalling(String referenceFile, List<String> bamFileList) {
		this.referenceFile = referenceFile;
		this.bamFileList = bamFileList;
	}

	public SVRecord run(SVRecord svRecord, Map<String, Integer> sampleIndex, QueryInterval[] queryIntervals)
			throws IOException {
		List<Integer> endList = new ArrayList();
		List<Integer> startList = new ArrayList();
		for (String bamFile : bamFileList) {
			File samReaderFile = new File(bamFile);
			try (SamReader samReader = SamReaderFactory.makeDefault().open(samReaderFile)) {
				String sample = samReader.getFileHeader().getReadGroups().get(0).getSample();
				SVType type = svRecord.getSVTypeArray()[sampleIndex.get(sample).intValue()];
				if (type == SVType.REFERENCE) {
					continue;
				}
				QueryInterval[] optimizedQueryIntervals = QueryInterval.optimizeIntervals(queryIntervals);
				SAMRecordIterator samRecordIterator = samReader.query(optimizedQueryIntervals, false);
				while (samRecordIterator.hasNext()) {
					SAMRecord samRecord = samRecordIterator.next();
					if (type == SVType.DELETION) {
						if (isLeftSoftClippedRead(samRecord) && samRecord.getReadPairedFlag()
								&& !samRecord.getMateUnmappedFlag()) {
							QueryInterval interval = optimizedQueryIntervals[optimizedQueryIntervals.length - 1];
							if (samRecord.getAlignmentStart() >= interval.start
									&& samRecord.getAlignmentStart() <= interval.end) {
								endList.add(samRecord.getAlignmentStart());
							}
						}
						if (isRightSoftClippedRead(samRecord) && samRecord.getReadPairedFlag()
								&& !samRecord.getMateUnmappedFlag()) {
							QueryInterval interval = optimizedQueryIntervals[0];
							if (samRecord.getAlignmentEnd() >= interval.start
									&& samRecord.getAlignmentEnd() <= interval.end) {
								startList.add(samRecord.getAlignmentEnd());
							}
						}
					}
					if (type == SVType.DUPLICATION) {
						if (isLeftSoftClippedRead(samRecord) && samRecord.getReadPairedFlag()
								&& !samRecord.getMateUnmappedFlag()) {
							QueryInterval interval = optimizedQueryIntervals[0];
							if (samRecord.getAlignmentStart() >= interval.start
									&& samRecord.getAlignmentStart() <= interval.end) {
								startList.add(samRecord.getAlignmentStart());
							}
						}

						if (isRightSoftClippedRead(samRecord) && samRecord.getReadPairedFlag()
								&& !samRecord.getMateUnmappedFlag()) {
							QueryInterval interval = optimizedQueryIntervals[optimizedQueryIntervals.length - 1];
							if (samRecord.getAlignmentEnd() >= interval.start
									&& samRecord.getAlignmentEnd() <= interval.end) {
								endList.add(samRecord.getAlignmentStart());
							}
						}
					}
				}
				samRecordIterator.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		int[] startRepeatElement = mostRepeatElement(startList);
		int[] endRepeatElement = mostRepeatElement(endList);
		if (startRepeatElement[1] > 1 && endRepeatElement[1] > 1) {
			int start = startRepeatElement[0];
			int end = endRepeatElement[0];
			if (start < end) {
				SVRecord svRecord2 = new SVRecord(svRecord.getChrom(), start, end, svRecord.getSVTypeArray(),
						svRecord.getEvidence() + ":BP");
				return svRecord2;
			} else {
				return null;
			}
		} else {
			return null;
		}
	}

	private int[] mostRepeatElement(List<Integer> elementList) {
		int[] arrayDemo = new int[elementList.size()];
		for (int i = 0; i < elementList.size(); i++) {
			arrayDemo[i] = elementList.get(i).intValue();
		}
		int[] element = new int[2];
		int result = -1;
		int N = arrayDemo.length;
		int majorityCount = 0;
		for (int i = 0; i < N; ++i) {
			int V = arrayDemo[i];
			int count = 0;
			for (int j = 0; j < N; ++j) {
				if (arrayDemo[j] == V) {
					count++;
				}
			}
			if (count > majorityCount) {
				majorityCount = count;
				result = V;
				element[0] = result;
				element[1] = majorityCount;
			}
		}
		return element;
	}

	private boolean isLeftSoftClippedRead(SAMRecord record) {
		if (record.getReadUnmappedFlag() || record.getCigar() == null || record.getMappingQuality() == 0) {
			return false;
		}
		List<CigarElement> elements = record.getCigar().getCigarElements();
		if (elements.get(0).getOperator() == CigarOperator.SOFT_CLIP) {
			return true;
		} else {
			return false;
		}
	}

	private boolean isRightSoftClippedRead(SAMRecord record) {
		if (record.getReadUnmappedFlag() || record.getCigar() == null || record.getMappingQuality() == 0) {
			return false;
		}
		List<CigarElement> elements = record.getCigar().getCigarElements();
		if (elements.get(elements.size() - 1).getOperator() == CigarOperator.SOFT_CLIP) {
			return true;
		} else {
			return false;
		}
	}
}
