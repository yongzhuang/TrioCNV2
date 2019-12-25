package cn.edu.hit.triocnv.discordantreadpair;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import cn.edu.hit.triocnv.util.Interval;
import cn.edu.hit.triocnv.util.IntervalReader;

/**
*
* @author Yongzhuang Liu
*/
public class ExcludingDRPs {

	private String drpFile;
	private String excludedRegion;

	public ExcludingDRPs(String drpFile, String excludedRegion) {
		super();
		this.drpFile = drpFile;
		this.excludedRegion = excludedRegion;
	}

	public List<ReadPair> exclude() throws IOException {
		List<ReadPair> readPairList = getReadPairList();
		if(excludedRegion==null || readPairList.size()==0) {
			return readPairList;
		}
		List<ReadPair> resultList = new ArrayList();
		String chrom = readPairList.get(0).getChrom();
		List<Interval> excludedIntervalList = getIntervalList(chrom);
		int i = 0;
		for (ReadPair readPair : readPairList) {
			boolean tag = true;
			for (int j = i; j < excludedIntervalList.size(); j++) {
				Interval currentInterval = excludedIntervalList.get(j);
				if (currentInterval.isOverlap(readPair.getFirstInterval())
						|| currentInterval.isOverlap(readPair.getSecondInterval())) {
					tag = false;
				}
				if (readPair.getLeftFirst() > currentInterval.getEnd()) {
					i = j;
				}
				if (readPair.getRightSecond() < currentInterval.getStart()) {
					break;
				}
			}
			if (tag) {
				resultList.add(readPair);
			}
		}
		return resultList;
	}

	public List<Interval> getIntervalList(String chrom) throws IOException {
		IntervalReader intervalReader = new IntervalReader(excludedRegion);
		List<Interval> intervalList = intervalReader.getIntervalList();
		List<Interval> resultList = new ArrayList();
		for (Interval interval : intervalList) {
			if (interval.getChrom().equals(chrom)) {
				resultList.add(interval);
			}
		}
		return resultList;
	}

	private List<ReadPair> getReadPairList() throws IOException {
		List<ReadPair> readPairList = new ArrayList();
		BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(drpFile)));
		String line = null;
		while ((line = bufferedReader.readLine()) != null) {
			String[] record = line.split("\t");
			String sample = record[0];
			String chrom = record[1];
			int leftFirst = Integer.parseInt(record[2]);
			int rightFirst = Integer.parseInt(record[3]);
			int leftSecond = Integer.parseInt(record[4]);
			int rightSecond = Integer.parseInt(record[5]);
			readPairList.add(new ReadPair(sample, chrom, leftFirst, rightFirst, leftSecond, rightSecond));
		}
		bufferedReader.close();
		return readPairList;
	}
}
