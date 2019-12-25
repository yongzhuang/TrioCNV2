package cn.edu.hit.triocnv.util;


/**
 *
 * @author Yongzhuang Liu
 */
public class Interval {

	private String chrom;
	private int start;
	private int end;

	public Interval(String chrom, int start, int end) {
		super();
		this.chrom = chrom;
		this.start = start;
		this.end = end;
	}

	public String getChrom() {
		return chrom;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}
	
	public boolean isOverlap(Interval interval) {
		if (!this.chrom.equals(interval.getChrom())) {
			return false;
		} else {
			if (this.getStart() > interval.getEnd() || this.getEnd() < interval.getStart()) {
				return false;
			} else
				return true;
		}
	}

	public String toString() {
		return chrom + ":" + start + "-" + end;
	}
}
