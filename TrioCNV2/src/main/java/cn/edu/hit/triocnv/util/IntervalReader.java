package cn.edu.hit.triocnv.util;

/**
*
* @author Yongzhuang Liu
*/

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class IntervalReader {

	String filename;

	public IntervalReader(String filename) {
		this.filename = filename;
	}

	public List<Interval> getIntervalList() throws IOException {
		List<Interval> intervalList = new ArrayList();
		BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(filename)));
		String line = null;
		while ((line = bufferedReader.readLine()) != null && line.trim().length() > 0) {
			String[] record = line.split("\t");
			String chrom = record[0];
			int start = Integer.parseInt(record[1]);
			int end = Integer.parseInt(record[2]);
			intervalList.add(new Interval(chrom, start, end));
		}
		bufferedReader.close();
		return intervalList;
	}
}
