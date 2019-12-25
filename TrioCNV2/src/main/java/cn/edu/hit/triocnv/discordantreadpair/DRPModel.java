package cn.edu.hit.triocnv.discordantreadpair;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.distribution.NormalDistribution;
import cn.edu.hit.triocnv.util.SVRecord;
import cn.edu.hit.triocnv.util.SVType;

/**
 *
 * @author Yongzhuang Liu
 */

public class DRPModel {

	private Set<ReadPair> readPairSet;
	private String fatherID;
	private String motherID;
	private String offspringID;
	private Map<String, int[]> parameterMap;

	public DRPModel(Set<ReadPair> readPairSet, Map<String, int[]> parameterMap, String fatherID, String motherID,
			String offspringID) {
		this.readPairSet = readPairSet;
		this.fatherID = fatherID;
		this.motherID = motherID;
		this.offspringID = offspringID;
		this.parameterMap = parameterMap;
	}

	public SVRecord call(SVType type) {
		int[] size = new int[] { getReadPairListOfSample(fatherID).size(), getReadPairListOfSample(motherID).size(),
				getReadPairListOfSample(offspringID).size() };
		if (size[2] > 1 && (size[0] > 1 || size[1] > 1)) {
			if (type == SVType.DELETION) {
				return getDEL();
			}
			if (type == SVType.DUPLICATION) {
				return getTandemDup();
			}
		}
		if (size[2] > 2 && size[0] < 2 && size[1] < 2) {
			if (type == SVType.DELETION) {
				return getDEL();
			}
			if (type == SVType.DUPLICATION) {
				return getTandemDup();
			}
		}
		if (size[2] < 2&& (size[0] > 2 || size[1] > 2)) {
			if (type == SVType.DELETION) {
				return getDEL();
			}
			if (type == SVType.DUPLICATION) {
				return getTandemDup();
			}
		}
		return null;
	}

	private List<ReadPair> getReadPairListOfSample(String sample) {
		List<ReadPair> readPairList = new ArrayList();
		for (ReadPair readPair : readPairSet) {
			if (readPair.getSample().equals(sample)) {
				readPairList.add(readPair);
			}
		}
		return readPairList;
	}

	public SVRecord getTandemDup() {
		int start = Integer.MAX_VALUE;
		int end = 0;
		ReadPair tmp = null;
		for (ReadPair readPair : readPairSet) {
			tmp = readPair;
			if (readPair.getRightFirst() < start) {
				start = readPair.getRightFirst();
			}
			if (readPair.getLeftSecond() > end) {
				end = readPair.getLeftSecond();
			}
		}
		int[] size = new int[] { getReadPairListOfSample(fatherID).size(), getReadPairListOfSample(motherID).size(),
				getReadPairListOfSample(offspringID).size() };

		SVType[] svTypeArray = new SVType[3];
		double quality = size[0] + size[1] + size[2];
		if (size[0] > 1) {
			svTypeArray[0] = SVType.DUPLICATION;
		} else {
			svTypeArray[0] = SVType.REFERENCE;
		}
		if (size[1] > 1) {
			svTypeArray[1] = SVType.DUPLICATION;
		} else {
			svTypeArray[1] = SVType.REFERENCE;
		}
		if (size[2] > 1) {
			svTypeArray[2] = SVType.DUPLICATION;
		} else {
			svTypeArray[2] = SVType.REFERENCE;
		}
		return new SVRecord(tmp.getChrom(), start, end, svTypeArray, quality);
	}

	public SVRecord getDEL() {
		int start = 0;
		int end = Integer.MAX_VALUE;
		ReadPair tmp = null;
		for (ReadPair readPair : readPairSet) {
			tmp = readPair;
			if (readPair.getRightFirst() > start) {
				start = readPair.getRightFirst();
			}
			if (readPair.getLeftSecond() < end) {
				end = readPair.getLeftSecond();
			}
		}
		int[] size = new int[] { getReadPairListOfSample(fatherID).size(), getReadPairListOfSample(motherID).size(),
				getReadPairListOfSample(offspringID).size() };
		double quality = getLikelihoodRatio();
		SVType[] svTypeArray = new SVType[3];
		if (size[0] > 1) {
			svTypeArray[0] = SVType.DELETION;
		} else {
			svTypeArray[0] = SVType.REFERENCE;
		}
		if (size[1] > 1) {
			svTypeArray[1] = SVType.DELETION;
		} else {
			svTypeArray[1] = SVType.REFERENCE;
		}
		if (size[2] > 1) {
			svTypeArray[2] = SVType.DELETION;
		} else {
			svTypeArray[2] = SVType.REFERENCE;
		}
		return new SVRecord(tmp.getChrom(), start, end, svTypeArray, quality);
	}

	private double getLikelihoodRatio() {

		int[] parameters = getMaxValuesOfParameters(parameterMap);
		int median = parameters[0];
		int sd = parameters[1];
		double prob1 = 0;
		NormalDistribution normalDistribution = new NormalDistribution(median, sd);
		for (ReadPair readPair : readPairSet) {
			int insertSize = readPair.getRightSecond() - readPair.getLeftFirst() + 1;
			prob1 += normalDistribution.logDensity(insertSize);
		}
		int start = 0;
		int end = Integer.MAX_VALUE;
		for (ReadPair readPair : readPairSet) {
			if (readPair.getRightFirst() > start) {
				start = readPair.getRightFirst();
			}
			if (readPair.getLeftSecond() < end) {
				end = readPair.getLeftSecond();
			}
		}
		int size = end - start + 1;
		double prob2 = 0;
		NormalDistribution delNormalDistribution = new NormalDistribution(median + size, sd);
		for (ReadPair readPair : readPairSet) {
			int insertSize = readPair.getRightSecond() - readPair.getLeftFirst() + 1;
			prob2 += delNormalDistribution.logDensity(insertSize);
		}
		double prob = prob2 - prob1;
		return prob;
	}

	private int[] getMaxValuesOfParameters(Map<String, int[]> parameterMap) {
		int[] parameters = new int[2];
		Set<String> sampleSet = parameterMap.keySet();
		int maxMedian = 0;
		int maxSD = 0;
		for (String sample : sampleSet) {
			int tmpMedian = parameterMap.get(sample)[0];
			int tmpSD = parameterMap.get(sample)[1];
			if (tmpMedian > maxMedian) {
				maxMedian = tmpMedian;
			}
			if (tmpSD > maxSD) {
				maxSD = tmpSD;
			}
		}
		parameters[0] = maxMedian;
		parameters[1] = maxSD;
		return parameters;
	}
}
