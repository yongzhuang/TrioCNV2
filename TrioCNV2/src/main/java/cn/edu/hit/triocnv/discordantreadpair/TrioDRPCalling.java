package cn.edu.hit.triocnv.discordantreadpair;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import cn.edu.hit.triocnv.util.CNVRecord;
import cn.edu.hit.triocnv.util.PEDReader;
import cn.edu.hit.triocnv.util.ParameterReader;
import cn.edu.hit.triocnv.util.SVRecord;
import cn.edu.hit.triocnv.util.SVType;
import cn.edu.hit.triocnv.util.Trio;

/**
 *
 * @author Yongzhuang Liu
 */

public class TrioDRPCalling {

	private String drpFolder;
	private String pedFile;
	private String outputFile;
	private String parameterFile;
	private String excludeFile;

	public TrioDRPCalling(String drpFolder, String pedFile, String parameterFile, String excludeFile,
			String outputFile) {
		super();
		this.drpFolder = drpFolder;
		this.pedFile = pedFile;
		this.outputFile = outputFile;
		this.parameterFile = parameterFile;
		this.excludeFile = excludeFile;
	}

	public void run() throws IOException {

		Trio trio = (new PEDReader(pedFile).getTrios()).get(0);
		String fatherID = trio.getFather().getIndividualID();
		String motherID = trio.getMother().getIndividualID();
		String offspringID = trio.getOffspring().getIndividualID();
		File drpDir = new File(drpFolder);
		ParameterReader parameterReader = new ParameterReader(parameterFile);
		Map<String, int[]> parameterMap = parameterReader.getParameters();
		List<SVRecord> svRecordList = new ArrayList();
		if (drpDir.exists() && drpDir.isDirectory()) {
			File[] files = drpDir.listFiles();
			for (File file : files) {
				String fileName = file.getAbsolutePath();
				ExcludingDRPs excluding = new ExcludingDRPs(fileName, excludeFile);
				List<ReadPair> readPairList = excluding.exclude();
				if (fileName.endsWith(".DEL")) {
					Clustering clustering = new Clustering(parameterMap, readPairList, fatherID, motherID, offspringID);
					//List<SVRecord> resultSVRecordList = merge(clustering.cluster(SVType.DELETION));
					List<SVRecord> resultSVRecordList = clustering.cluster(SVType.DELETION);
					if (resultSVRecordList != null) {
						svRecordList.addAll(resultSVRecordList);
					}
				}
				if (fileName.endsWith(".DUP")) {
					Clustering clustering = new Clustering(parameterMap, readPairList, fatherID, motherID, offspringID);
					//List<SVRecord> resultSVRecordList=merge(clustering.cluster(SVType.DUPLICATION));
					List<SVRecord> resultSVRecordList=clustering.cluster(SVType.DUPLICATION);
					if (resultSVRecordList != null) {
						svRecordList.addAll(resultSVRecordList);
					}
				}
			}
		}
		Collections.sort(svRecordList);
		FileWriter svWriter = new FileWriter(outputFile);
		svWriter.write("#CHROM" + "\t" + "START" + "\t" + "END" + "\t" + fatherID + "\t" + motherID + "\t" + offspringID
				+ "\t" + "Quality" + "\t" + "Evidence" + "\n");
		for (SVRecord svRecord : svRecordList) {
			if (svRecord.getLength() > 100) {
				svWriter.write(svRecord + "\tRP" + "\n");
			}
		}
		svWriter.close();
	}

	private List<SVRecord> merge(List<SVRecord> records) {
		if(records==null || records.size()==0) {
			return null;
		}
		List<SVRecord> result = new ArrayList<SVRecord>();
		List<SVRecord> cluster = new ArrayList<SVRecord>();
		cluster.add(records.get(0));
		for (int i = 1; i < records.size(); i++) {
			SVRecord svRecord = records.get(i);
			if (isOverlap(svRecord, cluster)) {
				cluster.add(svRecord);
			} else {
				result.add(getMaxQualitySVRecord(cluster));
				cluster = new ArrayList<SVRecord>();
				cluster.add(svRecord);
			}
		}
		result.add(getMaxQualitySVRecord(records));
		return result;
	}

	private boolean isOverlap(SVRecord svRecord, List<SVRecord> cluster) {
		for (SVRecord tmpSVRecord : cluster) {
			if (svRecord.isOverlap(tmpSVRecord)) {
				return true;
			}
		}
		return false;
	}

	private SVRecord getMaxQualitySVRecord(List<SVRecord> records) {
		if (records.size() == 0) {
			return null;
		}
		if (records.size() == 1) {
			return records.get(0);
		}
		double quality = 0;
		SVRecord svRecord = null;
		for (SVRecord tmpSVRecord : records) {
			if (tmpSVRecord.getQuality() > quality) {
				quality = tmpSVRecord.getQuality();
				svRecord = tmpSVRecord;
			}
		}
		return svRecord;
	}
}
