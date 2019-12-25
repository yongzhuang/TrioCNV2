package cn.edu.hit.triocnv.breakpoint;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.nio.channels.OverlappingFileLockException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import cn.edu.hit.triocnv.util.CNVRecord;
import cn.edu.hit.triocnv.util.PEDReader;
import cn.edu.hit.triocnv.util.ParameterReader;
import cn.edu.hit.triocnv.util.SVRecord;
import cn.edu.hit.triocnv.util.SVType;
import cn.edu.hit.triocnv.util.Trio;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqEncoder;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;

/**
 * 
 */

/**
 * @author Yongzhuang Liu
 *
 */
public class SingleThreadRefinement implements Runnable {
	private String referenceFile;
	private List<String> bamFileList;
	private String parameterFile;
	private List<SVRecord> svRecordList;
	private String outputFolder;
	private String pedFile;
	private int size;
	private int deviation;
	private static final int MAX_NUM_READS_ASSEMBLY=1000;

	public SingleThreadRefinement(List<SVRecord> svRecordList, String referenceFile, List<String> bamFileList,
			String pedFile, String parameterFile, String outputFolder, int size, int deviation) {
		this.svRecordList = svRecordList;
		this.bamFileList = bamFileList;
		this.referenceFile = referenceFile;
		this.outputFolder = outputFolder;
		this.parameterFile = parameterFile;
		this.pedFile = pedFile;
		this.size = size;
		this.deviation = deviation;
	}

	@Override
	public void run() {
		try {
			ParameterReader parameterReader = new ParameterReader(parameterFile);
			Map<String, int[]> parameterMap = parameterReader.getParameters();
			Trio trio = (new PEDReader(pedFile).getTrios()).get(0);
			String fatherID = trio.getFather().getIndividualID();
			String motherID = trio.getMother().getIndividualID();
			String offspringID = trio.getOffspring().getIndividualID();
			Map<String, Integer> sampleIndex = new HashMap();
			sampleIndex.put(fatherID, 0);
			sampleIndex.put(motherID, 1);
			sampleIndex.put(offspringID, 2);
			File outputDir = new File(outputFolder);
			if (outputDir.exists()) {
				if (!outputDir.isDirectory()) {
					return;
				}
			} else {
				outputDir.mkdir();
			}
			List<SVRecord> resultSVRecordList = new ArrayList();
			for (int i = 0; i < svRecordList.size(); i++) {
				System.out.println("i=" + i);
				SVRecord svRecord = svRecordList.get(i);
				int referenceIndex = -1;
				String chrom = svRecord.getChrom();
				int start = svRecord.getStart();
				int end = svRecord.getEnd();
				if (chrom.equals("chrX") || chrom.equals("X")) {
					referenceIndex = 22;
				} else if (chrom.equals("chrY") || chrom.equals("Y")) {
					referenceIndex = 23;
				} else {
					if (chrom.startsWith("chr")) {
						referenceIndex = Integer.parseInt(chrom.substring(3)) - 1;
					} else {
						referenceIndex = Integer.parseInt(chrom) - 1;
					}
				}
				QueryInterval[] queryIntervals = null;
				int inner = size / 2;
				int outer = size / 2;
				if ((end - start + 1) > 2 * inner) {
					queryIntervals = new QueryInterval[2];
					QueryInterval leftQueryInterval = new QueryInterval(referenceIndex, start - outer, start + inner);
					QueryInterval rightQueryInterval = new QueryInterval(referenceIndex, end - inner, end + outer);
					queryIntervals[0] = leftQueryInterval;
					queryIntervals[1] = rightQueryInterval;
				} else {
					queryIntervals = new QueryInterval[1];
					QueryInterval queryInterval = new QueryInterval(referenceIndex, start - outer, end + outer);
					queryIntervals[0] = queryInterval;
				}
				queryIntervals = QueryInterval.optimizeIntervals(queryIntervals);
				SVRecord softclipSVRecord = new SoftClipCalling(referenceFile, bamFileList).run(svRecord, sampleIndex,
						queryIntervals);
				if (softclipSVRecord != null) {
					resultSVRecordList.add(softclipSVRecord);
				} else {
					String prefixFile = svRecord.getChrom() + ":" + svRecord.getStart() + "-" + svRecord.getEnd();
					String assemblyWorkDir = outputFolder + "/" + prefixFile;
					File assemblyWorkDirFile = new File(assemblyWorkDir);
					if (assemblyWorkDirFile.exists()) {
						if (!assemblyWorkDirFile.isDirectory()) {
							return;
						}
					} else {
						assemblyWorkDirFile.mkdir();
					}
					String fastqFile1 = assemblyWorkDirFile.getAbsolutePath() + "/" + prefixFile + "_1.fastq";
					String fastqFile2 = assemblyWorkDirFile.getAbsolutePath() + "/" + prefixFile + "_2.fastq";
					extractSVReads(fastqFile1, fastqFile2, svRecord, parameterMap, sampleIndex, queryIntervals,
							deviation);
					int numOfReads = 0;
					if ((new File(fastqFile1)).exists()) {
						FastqReader fastqReader = new FastqReader(new File(fastqFile1));
						Iterator<FastqRecord> iterator = fastqReader.iterator();
						while (iterator.hasNext()) {
							iterator.next();
							numOfReads++;
						}
						fastqReader.close();
					}
					if (numOfReads > 1 && numOfReads < MAX_NUM_READS_ASSEMBLY) {
						Assembly assembly = new Assembly(fastqFile1, fastqFile2, assemblyWorkDir, 1);
						assembly.runCommand();
						String contigFile = assemblyWorkDirFile.getAbsolutePath() + "/" + "primary-contigs.fa";
						if ((new File(contigFile)).exists()) {
							BufferedReader reader = new BufferedReader(new FileReader(new File(contigFile)));
							String line = reader.readLine();
							reader.close();
							if (line != null) {
								String fastaFile = assemblyWorkDirFile.getAbsolutePath() + "/" + "reference.fa";
								ExtractSequence extractSequence = new ExtractSequence(referenceFile, fastaFile,
										svRecord.getChrom());
								extractSequence.getSubSequence(chrom, svRecord.getStart() - outer - 1000,
										svRecord.getEnd() + outer + 1000);
								Alignment alignment = new Alignment(contigFile, fastaFile,
										assemblyWorkDirFile.getAbsolutePath());
								alignment.alignContig();
								String samFile = assemblyWorkDirFile.getAbsolutePath() + "/" + "primary-contigs.sam";
								AlignmentParser alignmentParser = new AlignmentParser(svRecord, samFile, outer);
								SVRecord assemblySVRecord = null;
								SVType[] type = svRecord.getSVTypeArray();
								if (type[0] == SVType.DELETION || type[1] == SVType.DELETION
										|| type[2] == SVType.DELETION) {
									assemblySVRecord = alignmentParser.parse();
								} else {
									assemblySVRecord = alignmentParser.parseTandemDuplication();
								}
								if (assemblySVRecord != null) {
									resultSVRecordList.add(assemblySVRecord);
								}
							}
						}
					}
					deleteDir(assemblyWorkDirFile);
				}
			}
			FileOutputStream out = new FileOutputStream(new File(outputFolder + "/Temp_CNV.txt"), true);
			FileChannel fcout = out.getChannel();
			FileLock flout = null;
			while (true) {
				try {
					flout = fcout.tryLock();
					break;
				} catch (OverlappingFileLockException e) {
					try {
						Thread.sleep(100);
					} catch (InterruptedException ex) {
						ex.printStackTrace();
					}
				}
			}
			for (SVRecord tmpSVRecord : resultSVRecordList) {
				String tmp = tmpSVRecord.getChrom() + "\t" + tmpSVRecord.getStart() + "\t" + tmpSVRecord.getEnd()
						+ "\t" + tmpSVRecord.getSVTypeArray()[0] + "\t" + tmpSVRecord.getSVTypeArray()[1] + "\t"
						+ tmpSVRecord.getSVTypeArray()[2] + "\t" + tmpSVRecord.getEvidence() + "\n";
				out.write(tmp.getBytes());
			}
			flout.release();
			fcout.close();
			out.close();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	private void extractSVReads(String fastqFile1, String fastqFile2, SVRecord svRecord,
			Map<String, int[]> parameterMap, Map<String, Integer> sampleIndex, QueryInterval[] queryIntervals,
			int deviation) throws IOException {
		FastqWriter fastqWriter1 = null;
		FastqWriter fastqWriter2 = null;
		for (String bamFile : bamFileList) {
			SamReader samReader = SamReaderFactory.makeDefault().open(new File(bamFile));
			String sample = samReader.getFileHeader().getReadGroups().get(0).getSample();
			SVType type = svRecord.getSVTypeArray()[sampleIndex.get(sample).intValue()];
			if (type == SVType.REFERENCE) {
				continue;
			}
			int[] parameters = parameterMap.get(sample);
			int median = parameters[0];
			int sd = parameters[1];
			HashMap<String, SAMRecord> samRecordTable = new HashMap<String, SAMRecord>();
			SAMRecordIterator samRecordIterator = samReader.query(queryIntervals, true);
			int numOfReads = 0;
			while (samRecordIterator.hasNext()) {
				SAMRecord samRecord = samRecordIterator.next();
				String readName = samRecord.getReadName();
				if (samRecord.getReadPairedFlag() && !samRecord.isSecondaryOrSupplementary()) {
					if (samRecordTable.containsKey(readName)) {
						SAMRecord samRecord2 = samRecordTable.get(readName);
						if ((type == SVType.DELETION && isDeletionDRPs(samRecord, samRecord2, median, sd, deviation))
								|| (type == SVType.DUPLICATION && isDuplicationDRPs(samRecord, samRecord2))
								|| isOneEndAnchoredRead(samRecord, samRecord2)
								|| isSoftClippedRead(samRecord, samRecord2)) {
							if (fastqWriter1 == null) {
								fastqWriter1 = new BasicFastqWriter(new File(fastqFile1));
							}
							if (fastqWriter2 == null) {
								fastqWriter2 = new BasicFastqWriter(new File(fastqFile2));
							}
							FastqRecord fastqRecord = FastqEncoder.asFastqRecord(samRecord);
							FastqRecord fastqRecord2 = FastqEncoder.asFastqRecord(samRecord2);
							if (isProperPairedReads(samRecord, samRecord2)) {
								numOfReads++;
								if(numOfReads>MAX_NUM_READS_ASSEMBLY) {
									if (fastqWriter1 != null && fastqWriter2 != null) {
										fastqWriter1.close();
										fastqWriter2.close();
									}
									return;
								}
								if (samRecord.getFirstOfPairFlag()) {
									fastqWriter1.write(fastqRecord);
								}
								if (samRecord.getSecondOfPairFlag()) {
									fastqWriter2.write(fastqRecord);
								}
								if (samRecord2.getFirstOfPairFlag()) {
									fastqWriter1.write(fastqRecord2);
								}
								if (samRecord2.getSecondOfPairFlag()) {
									fastqWriter2.write(fastqRecord2);
								}
							}
						}
						samRecordTable.remove(readName);
					} else {
						samRecordTable.put(readName, samRecord);
					}
				}
			}
			samRecordIterator.close();
			samReader.close();

			if (samRecordTable.size() > 0) {
				List<QueryInterval> mateIntervalList = new ArrayList();
				Set<String> keySet = samRecordTable.keySet();
				Iterator<String> iter = keySet.iterator();
				while (iter.hasNext()) {
					String key = iter.next();
					SAMRecord samRecord = samRecordTable.get(key);
					if (samRecord.getMateReferenceIndex() != -1) {
						mateIntervalList.add(new QueryInterval(samRecord.getMateReferenceIndex(),
								samRecord.getMateAlignmentStart(), samRecord.getMateAlignmentStart()));
					}
				}
				Collections.sort(mateIntervalList);
				QueryInterval[] mateIntervals = new QueryInterval[mateIntervalList.size()];
				mateIntervalList.toArray(mateIntervals);
				HashMap<String, SAMRecord> mateSAMRecordTable = getMateSAMRecordTable(mateIntervals, samRecordTable,
						bamFile);
				Set<String> keySet2 = samRecordTable.keySet();
				Iterator<String> iter2 = keySet2.iterator();
				while (iter2.hasNext()) {
					String key = iter2.next();
					SAMRecord samRecord = samRecordTable.get(key);
					SAMRecord mateSAMRecord = mateSAMRecordTable.get(key);
					if (samRecord == null || mateSAMRecord == null) {
						continue;
					}
					if ((type == SVType.DELETION && isDeletionDRPs(samRecord, mateSAMRecord, median, sd, deviation))
							|| (type == SVType.DUPLICATION && isDuplicationDRPs(samRecord, mateSAMRecord))
							|| isOneEndAnchoredRead(samRecord, mateSAMRecord)
							|| isSoftClippedRead(samRecord, mateSAMRecord)) {
						FastqRecord fastqRecord = FastqEncoder.asFastqRecord(samRecord);
						FastqRecord mateFastqRecord = FastqEncoder.asFastqRecord(mateSAMRecord);

						if (fastqWriter1 == null) {
							fastqWriter1 = new BasicFastqWriter(new File(fastqFile1));
						}
						if (fastqWriter2 == null) {
							fastqWriter2 = new BasicFastqWriter(new File(fastqFile2));
						}
						if (isProperPairedReads(samRecord, mateSAMRecord)) {
							numOfReads++;
							if (samRecord.getFirstOfPairFlag()) {
								fastqWriter1.write(fastqRecord);
							}
							if (samRecord.getSecondOfPairFlag()) {
								fastqWriter2.write(fastqRecord);
							}
							if (mateSAMRecord.getFirstOfPairFlag()) {
								fastqWriter1.write(mateFastqRecord);
							}
							if (mateSAMRecord.getSecondOfPairFlag()) {
								fastqWriter2.write(mateFastqRecord);
							}
							if (numOfReads > MAX_NUM_READS_ASSEMBLY) {
								if (fastqWriter1 != null && fastqWriter2 != null) {
									fastqWriter1.close();
									fastqWriter2.close();
								}
								return;
							}
						}
					}
				}
			}
		}
		if (fastqWriter1 != null && fastqWriter2 != null) {
			fastqWriter1.close();
			fastqWriter2.close();
		}
	}

	private HashMap<String, SAMRecord> getMateSAMRecordTable(QueryInterval[] queryIntervals,
			HashMap<String, SAMRecord> samRecordTable, String bamFile) throws IOException {
		QueryInterval[] intervals = QueryInterval.optimizeIntervals(queryIntervals);
		SamReader samReader = SamReaderFactory.makeDefault().open(new File(bamFile));
		SAMRecordIterator samRecordIterator = samReader.query(intervals, false);
		HashMap<String, SAMRecord> mateRecordTable = new HashMap<String, SAMRecord>();
		while (samRecordIterator.hasNext()) {
			SAMRecord mateRecord = samRecordIterator.next();
			String readName = mateRecord.getReadName();
			if (samRecordTable.containsKey(readName)
					&& !samRecordTable.get(readName).getReadString().equals(mateRecord.getReadString())) {
				mateRecordTable.put(mateRecord.getReadName(), mateRecord);
			}
		}
		samRecordIterator.close();
		samReader.close();
		return mateRecordTable;
	}

	private boolean isOneEndAnchoredRead(SAMRecord record, SAMRecord record2) {
		boolean oeaTag = record.getReadPairedFlag() && (record.getReadUnmappedFlag() ^ record.getMateUnmappedFlag());
		boolean oeaTag2 = record2.getReadPairedFlag()
				&& (record2.getReadUnmappedFlag() ^ record2.getMateUnmappedFlag());
		return oeaTag || oeaTag2;
	}

	private boolean isSplitRead(SAMRecord record, SAMRecord record2) {
		return record.getAttribute(SAMTag.SA.name()) != null || record2.getAttribute(SAMTag.SA.name()) != null;
	}

	private boolean isSoftClippedRead(SAMRecord record) {
		if (record.getReadUnmappedFlag() || record.getCigar() == null) {
			return false;
		}
		List<CigarElement> elements = record.getCigar().getCigarElements();
		int i = 0;
		while (i < elements.size()) {
			if (elements.get(i).getOperator() == CigarOperator.SOFT_CLIP) {
				return true;
			}
			i++;
		}
		return false;
	}

	private boolean isSoftClippedRead(SAMRecord record, SAMRecord record2) {

		if (isSoftClippedRead(record) || isSoftClippedRead(record2)) {
			return true;
		} else {
			return false;
		}
	}

	private boolean isDeletionDRPs(SAMRecord record, int median, int sd, int deviation) {
		if (!record.getReadPairedFlag() || record.getReadUnmappedFlag() || record.getMateUnmappedFlag()
				|| record.getReferenceIndex() != record.getMateReferenceIndex()) {
			return false;
		}
		SamPairUtil.PairOrientation orientation = SamPairUtil.getPairOrientation(record);
		if (orientation == SamPairUtil.PairOrientation.FR) {
			int insertSize = Math.abs(record.getInferredInsertSize());
			if (insertSize > median + sd * deviation) {
				return true;
			}
		}
		return false;
	}

	private boolean isDeletionDRPs(SAMRecord record, SAMRecord record2, int median, int sd, int deviation) {
		if (isDeletionDRPs(record, median, sd, deviation) && isDeletionDRPs(record2, median, sd, deviation)) {
			return true;
		} else {
			return false;
		}
	}

	private boolean isDuplicationDRPs(SAMRecord record) {
		if (!record.getReadPairedFlag() || record.getReadUnmappedFlag() || record.getMateUnmappedFlag()
				|| record.getReferenceIndex() != record.getMateReferenceIndex()) {
			return false;
		}
		SamPairUtil.PairOrientation orientation = SamPairUtil.getPairOrientation(record);
		if (orientation == SamPairUtil.PairOrientation.RF) {
			return true;
		}
		return false;
	}

	private boolean isDuplicationDRPs(SAMRecord record, SAMRecord record2) {
		if (isDuplicationDRPs(record) && isDuplicationDRPs(record2)) {
			return true;
		} else {
			return false;
		}
	}

	private boolean deleteDir(File dir) {
		if (!dir.exists())
			return false;
		if (dir.isDirectory()) {
			String[] childrens = dir.list();
			for (String child : childrens) {
				boolean success = deleteDir(new File(dir, child));
				if (!success)
					return false;
			}
		}
		return dir.delete();
	}

	private boolean isProperPairedReads(SAMRecord record1, SAMRecord record2) {
		if (record1.getReadName().equals(record2.getReadName())) {
			if ((record1.getFirstOfPairFlag() && record2.getSecondOfPairFlag())
					|| (record1.getSecondOfPairFlag() && record2.getFirstOfPairFlag())) {
				return true;
			} else {
				return false;
			}
		} else {
			return false;
		}
	}

}
