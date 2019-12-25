package cn.edu.hit.triocnv.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
*
* @author Yongzhuang Liu
*/

public class ParameterReader {
	
	String filename;

	public ParameterReader(String filename) {
		this.filename = filename;
	}

	public Map<String, int[]> getParameters() throws IOException {
		Map<String, int[]> parameterMap = new HashMap();
		BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(filename)));
		String line = null;
		while ((line = bufferedReader.readLine()) != null && line.trim().length() > 0) {
			String[] record = line.split("\t");
			int[] parameters = new int[2];
			String sample=record[0];
			parameters[0] =(int)Double.parseDouble(record[1]);
			parameters[1] = (int)Double.parseDouble(record[2]);
			parameterMap.put(sample, parameters);
		}
		bufferedReader.close();
		return parameterMap;
	}

}
