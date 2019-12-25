package cn.edu.hit.triocnv.readdepth;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.PascalDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;

import cn.edu.hit.triocnv.readdepth.EmissionException;

/**
 *
 * @author Yongzhuang Liu
 */
public class Emission {

	private NBModel[] nbModels;

	private static final double LOG_MIN_PRO = -1000;
	private static final double LOG_MAX_PRO = -Math.pow(10, -100);
	private static final int MAX_RD = 5000;
	private static final double EPSILON = 0.1;
	private static final double DISPERSION_THRESHOLD = 500;
	private static final int MIN_READ_COUNT = 1;

	public Emission(NBModel[] nbModels) {
		this.nbModels = nbModels;
	}

	public double getProbOfRDEmission(int[] states, Observation observation) throws Exception {
		double prob = 0;
		int[] rd = observation.getRD();
		int gc = observation.getGC();
		double mappability = observation.getMappability();
		for (int i = 0; i < 3; i++) {
			prob += getProbOfRDGivenState(states[i], gc, mappability, rd[i], nbModels[i]);
		}
		if (prob == Double.POSITIVE_INFINITY || prob == Double.NEGATIVE_INFINITY || Double.isNaN(prob)) {
			throw new EmissionException();
		} else {
			return prob;
		}
	}

	private double getProbOfRDGivenState(int state, int gc, double mappability, int rd, NBModel nbModel) {
		double[] theta = nbModel.getTheta();
		double[] coef1 = nbModel.getCoef1();
		double[] coef2 = nbModel.getCoef2();
		if (theta[gc] == 0) {
			return 0;
		}
		double prob = 0;
		if (rd == 0) {
			if (state == 0) {
				prob = LOG_MAX_PRO;
			} else {
				prob = LOG_MIN_PRO;
			}
		} else if (rd > MAX_RD) {
			if (state == 4) {
				prob = LOG_MAX_PRO;
			} else {
				prob = LOG_MIN_PRO;
			}
		} else {
			double mu = Math.exp(coef1[gc] + coef2[gc] * mappability);
			double size = theta[gc];
			double p = 0;
			if (state == 0) {
				size = size * EPSILON;
				if (size > DISPERSION_THRESHOLD) {
					size = DISPERSION_THRESHOLD;
				}
				p = size / (size + mu * EPSILON);
			} else {
				size = size * (0.5 * (double) state);
				if (size > DISPERSION_THRESHOLD) {
					size = DISPERSION_THRESHOLD;
				}
				p = size / (size + mu * 0.5 * (double) state);
			}
			PascalDistribution pascal = new PascalDistribution((int) Math.ceil(size), p);
			prob = pascal.logProbability(rd);
			if (prob == Double.NEGATIVE_INFINITY) {
				prob = LOG_MIN_PRO;
			}
		}
		return prob;
	}
}
