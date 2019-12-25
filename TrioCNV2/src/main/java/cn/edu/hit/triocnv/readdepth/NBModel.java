package cn.edu.hit.triocnv.readdepth;

/**
 *
 * @author Yongzhuang Liu
 */
public class NBModel {

	private double[] theta;
	private double[] coef1;
	private double[] coef2;

	public NBModel(double[] theta, double[] coef1, double[] coef2) {
		this.theta = theta;
		this.coef1 = coef1;
		this.coef2 = coef2;
	}

	public double[] getTheta() {
		return theta;
	}

	public void setTheta(double[] theta) {
		this.theta = theta;
	}

	public double[] getCoef1() {
		return coef1;
	}

	public void setCoef1(double[] coef1) {
		this.coef1 = coef1;
	}

	public double[] getCoef2() {
		return coef2;
	}

	public void setCoef2(double[] coef2) {
		this.coef2 = coef2;
	}
}
