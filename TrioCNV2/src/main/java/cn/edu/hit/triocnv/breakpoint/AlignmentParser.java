package cn.edu.hit.triocnv.breakpoint;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import cn.edu.hit.triocnv.util.SVRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMTag;

/**
 *
 * @author Yongzhuang Liu
 */

public class AlignmentParser {

	private String samFile;
	private SVRecord svRecord;
	private int outer;

	public AlignmentParser(SVRecord svRecord, String samFile, int outer) {
		this.svRecord = svRecord;
		this.samFile = samFile;
		this.outer = outer;
	}

	public SVRecord parse() throws IOException {
		int[] interval = new int[2];
		HashMap<String, SAMRecord> samRecordTable = new HashMap<String, SAMRecord>();

		SVRecord resultSVRecord = null;
		try (SamReader samReader = SamReaderFactory.makeDefault().open(new File(samFile))) {
			SAMRecordIterator samRecordIterator = samReader.iterator();
			while (samRecordIterator.hasNext()) {
				SAMRecord samRecord = samRecordIterator.next();
				if (samRecord.getMappingQuality() > 0) {
					String readName = samRecord.getReadName();
					if (samRecordTable.containsKey(readName)) {
						if (isRightClippedRead(samRecord)) {
							interval[0] = samRecord.getAlignmentEnd() + 1;
						}
						if (isLeftClippedRead(samRecord)) {
							interval[1] = samRecord.getAlignmentStart() - 1;
						}
						SAMRecord samRecord2 = samRecordTable.get(readName);
						if (isRightClippedRead(samRecord2)) {
							interval[0] = samRecord2.getAlignmentEnd() + 1;
						}
						if (isLeftClippedRead(samRecord2)) {
							interval[1] = samRecord2.getAlignmentStart() - 1;
						}
						samRecordTable.remove(readName);
						if (interval[0] > 0 && interval[1] > 0) {
							int start = svRecord.getStart() - outer - 1000 - 1 + interval[0];
							int end = svRecord.getStart() - outer - 1000 - 1 + interval[1];
							if (start < end) {
								resultSVRecord = new SVRecord(svRecord.getChrom(), start, end,
										svRecord.getSVTypeArray(), svRecord.getEvidence() + ":BP");
								break;
							}
						}
					} else {
						samRecordTable.put(readName, samRecord);
					}
				}
			}
			samRecordIterator.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return resultSVRecord;
	}

	public SVRecord parseTandemDuplication() throws IOException {
		int[] interval = new int[2];
		HashMap<String, SAMRecord> samRecordTable = new HashMap<String, SAMRecord>();
		SVRecord resultSVRecord = null;
		try (SamReader samReader = SamReaderFactory.makeDefault().open(new File(samFile))) {
			SAMRecordIterator samRecordIterator = samReader.iterator();
			while (samRecordIterator.hasNext()) {
				SAMRecord samRecord = samRecordIterator.next();
				if (samRecord.getMappingQuality() > 0) {
					String readName = samRecord.getReadName();
					if (samRecordTable.containsKey(readName)) {
						if (isLeftClippedRead(samRecord)) {
							interval[0] = samRecord.getAlignmentStart();
						}
						if (isRightClippedRead(samRecord)) {
							interval[1] = samRecord.getAlignmentEnd();
						}
						SAMRecord samRecord2 = samRecordTable.get(readName);
						if (isLeftClippedRead(samRecord2)) {
							interval[0] = samRecord2.getAlignmentStart();
						}
						if (isRightClippedRead(samRecord2)) {
							interval[1] = samRecord2.getAlignmentEnd();
						}
						if (interval[0] > 0 && interval[1] > 0) {
							int start = svRecord.getStart() - outer - 1000 - 1 + interval[0];
							int end = svRecord.getStart() - outer - 1000 - 1 + interval[1];
							if (start < end) {
								resultSVRecord = new SVRecord(svRecord.getChrom(), start, end,
										svRecord.getSVTypeArray(), svRecord.getEvidence() + ":BP");
								break;
							}
						}
						samRecordTable.remove(readName);
					} else {
						samRecordTable.put(readName, samRecord);
					}
				}
			}
			samRecordIterator.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return resultSVRecord;
	}

	private boolean isLeftClippedRead(SAMRecord record) {
		if (record.getReadUnmappedFlag() || record.getCigar() == null) {
			return false;
		}
		List<CigarElement> elements = record.getCigar().getCigarElements();
		if (elements.get(0).getOperator() == CigarOperator.S || elements.get(0).getOperator() == CigarOperator.H) {
			return true;
		} else {
			return false;
		}
	}

	private boolean isRightClippedRead(SAMRecord record) {
		if (record.getReadUnmappedFlag() || record.getCigar() == null) {
			return false;
		}
		List<CigarElement> elements = record.getCigar().getCigarElements();
		if (elements.get(elements.size() - 1).getOperator() == CigarOperator.S
				|| elements.get(elements.size() - 1).getOperator() == CigarOperator.H) {
			return true;
		} else {
			return false;
		}
	}
}
