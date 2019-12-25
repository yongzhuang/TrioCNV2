package cn.edu.hit.triocnv.discordantreadpair;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jgrapht.Graph;
import org.jgrapht.alg.clique.BronKerboschCliqueFinder;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import cn.edu.hit.triocnv.util.SVRecord;
import cn.edu.hit.triocnv.util.SVType;

/**
 *
 * @author Yongzhuang Liu
 */

public class Clustering {

	Map<String, int[]> parameterMap;
	private List<ReadPair> readPairList;
	private String fatherID;
	private String motherID;
	private String offspringID;

	private final static int MAX_READ_PAIRS = 500;
	private final static int MIN_READ_PAIRS = 2;
	private final static int DEVIATION = 3;

	public Clustering(Map<String, int[]> parameterMap, List<ReadPair> readPairList, String fatherID, String motherID,
			String offspringID) {
		this.parameterMap = parameterMap;
		this.readPairList = readPairList;
		this.fatherID = fatherID;
		this.motherID = motherID;
		this.offspringID = offspringID;
	}

	private SVRecord getClique(Graph<ReadPair, DefaultEdge> g, SVType type) {
		BronKerboschCliqueFinder clique = new BronKerboschCliqueFinder(g);
		Iterator<Set<ReadPair>> iterator = clique.maximumIterator();
		double maxQuality = -1;
		SVRecord svRecord = null;
		while (iterator.hasNext()) {
			Set<ReadPair> maxClique = (Set<ReadPair>) iterator.next();
			DRPModel trioDRPCalling = new DRPModel(maxClique, parameterMap, fatherID, motherID, offspringID);
			SVRecord tmpRecord = trioDRPCalling.call(type);
			if (tmpRecord != null && tmpRecord.getQuality() > maxQuality) {
				svRecord = tmpRecord;
			}
		}
		return svRecord;
	}

	public List<SVRecord> cluster(SVType type) throws IOException {
		int[] parameters = getMaxValuesOfParameters(parameterMap);
		int median = parameters[0];
		int standardDeviation = parameters[1];
		List<SVRecord> outputRecordList = new ArrayList();
		int[] visit = new int[readPairList.size()];
		List<ReadPair> nodeList = null;
		ReadPair lastReadPair = null;
		for (int i = 0; i < readPairList.size(); i++) {
			if (visit[i] == 0) {
				nodeList = new ArrayList();
				nodeList.add(readPairList.get(i));
				visit[i] = 1;
				lastReadPair = readPairList.get(i);
				for (int j = i + 1; j < readPairList.size(); j++) {
					if (visit[j] == 0) {
						ReadPair readPair = readPairList.get(j);
						if (Math.abs(readPair.getLeftFirst() - lastReadPair.getLeftFirst()) > (median
								+ DEVIATION * standardDeviation)) {
							break;
						}
						for (int k = nodeList.size() - 1; k >= 0; k--) {
							ReadPair node = nodeList.get(k);
							if (Math.abs(readPair.getLeftFirst() - node.getLeftFirst()) < (median
									+ DEVIATION * standardDeviation)
									&& Math.abs(readPair.getRightSecond() - node.getRightSecond()) < (median
											+ DEVIATION * standardDeviation)) {
								if (type == SVType.DELETION) {
									if (Math.max(readPair.getRightFirst(), node.getRightFirst()) < Math
											.min(readPair.getLeftSecond(), node.getLeftSecond())) {
										nodeList.add(readPair);
										visit[j] = 1;
										lastReadPair = readPair;
										break;
									}
								}
								if (type == SVType.DUPLICATION) {
									nodeList.add(readPair);
									visit[j] = 1;
									lastReadPair = readPair;
									break;
								}
							}
						}
					}
				}
				if (nodeList.size() < MAX_READ_PAIRS && nodeList.size() > MIN_READ_PAIRS) {
					Graph<ReadPair, DefaultEdge> g = addReadPairEages(type,nodeList, median, standardDeviation);
					SVRecord outputRecord = getClique(g, type);
					if (outputRecord != null) {
						outputRecordList.add(outputRecord);
					}
				}
			}
		}
		return outputRecordList;
	}

	public Graph<ReadPair, DefaultEdge> addReadPairEages(SVType type, List<ReadPair> readPairList, int median,
			int standardDeviation) {
		List<ReadPair> nodeList = new ArrayList();
		Graph<ReadPair, DefaultEdge> g = new SimpleGraph<>(DefaultEdge.class);
		for (int i = 0; i < readPairList.size(); i++) {
			if (!g.containsVertex(readPairList.get(i))) {
				g.addVertex(readPairList.get(i));
				nodeList.add(readPairList.get(i));
			}
		}
		for (int i = 0; i < nodeList.size(); i++) {
			ReadPair startNode = nodeList.get(i);
			for (int j = i + 1; j < nodeList.size(); j++) {
				ReadPair endNode = nodeList.get(j);
				if (Math.abs(
						endNode.getLeftFirst() - startNode.getLeftFirst()) < (median + DEVIATION * standardDeviation)
						&& Math.abs(endNode.getRightSecond() - startNode.getRightSecond()) < (median
								+ DEVIATION * standardDeviation)) {
					if (type == SVType.DELETION) {
						if (Math.max(endNode.getRightFirst(), endNode.getRightFirst()) < Math
								.min(endNode.getLeftSecond(), endNode.getLeftSecond())) {
							g.addEdge(startNode, endNode);
						}
					}
					if (type == SVType.DUPLICATION) {
						g.addEdge(startNode, endNode);
					}
				}
			}
		}
		return g;
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