package cn.edu.hit.triocnv.discordantreadpair;

import cn.edu.hit.triocnv.util.Interval;

/**
*
* @author Yongzhuang Liu
*/
public class ReadPair implements Comparable<ReadPair> {

	private String sample;
	private String chrom;
	private int leftFirst;
	private int rightFirst;
	private int leftSecond;
	private int rightSecond;

	public ReadPair(String sample, String chrom, int leftFirst, int rightFirst, int leftSecond, int rightSecond) {
		this.sample = sample;
		this.chrom = chrom;
		this.leftFirst = leftFirst;
		this.rightFirst = rightFirst;
		this.leftSecond = leftSecond;
		this.rightSecond = rightSecond;
	}

	public String getSample() {
		return sample;
	}

	public String getChrom() {
		return chrom;
	}

	public int getLeftFirst() {
		return leftFirst;
	}

	public int getRightFirst() {
		return rightFirst;
	}

	public int getLeftSecond() {
		return leftSecond;
	}

	public int getRightSecond() {
		return rightSecond;
	}

	public Interval getFirstInterval() {
		return new Interval(chrom, leftFirst, rightFirst);
	}

	public Interval getSecondInterval() {
		return new Interval(chrom, leftSecond, rightSecond);
	}

	@Override
	public String toString() {
		return this.sample + "\t" + this.chrom + "\t" + this.leftFirst + "\t" + this.rightFirst + "\t" + this.leftSecond
				+ "\t" + this.rightSecond;
	}

	@Override
	public int compareTo(ReadPair readPair) {
		if (!this.getChrom().equals(readPair.getChrom())) {
			return this.getChrom().compareTo(readPair.getChrom());
		} else {
			if (this.getLeftFirst() > readPair.getLeftFirst()) {
				return 1;
			} else if (this.getLeftFirst() < readPair.getLeftFirst()) {
				return -1;
			} else {
				if (this.getLeftSecond() > readPair.getLeftSecond()) {
					return 1;
				} else if (this.getLeftSecond() < readPair.getLeftSecond()) {
					return -1;
				} else {
					return 0;
				}
			}
		}
	}
}
