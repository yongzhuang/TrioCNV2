package cn.edu.hit.triocnv.util;

/**
*
* @author Yongzhuang Liu
*/
/**
 * @author Yongzhuang Liu
 *
 */
public class SVRecord implements Comparable<SVRecord> {

	private String chrom;
	private int start;
	private int end;
	private SVType[] svTypeArray;
	private double quality;
	private String evidence;

	public SVRecord(String chrom, int start, int end, SVType[] svTypeArray) {
		super();
		this.chrom = chrom;
		this.start = start;
		this.end = end;
		this.svTypeArray = svTypeArray;
	}

	public SVRecord(String chrom, int start, int end, SVType[] svTypeArray, double quality) {
		super();
		this.chrom = chrom;
		this.start = start;
		this.end = end;
		this.svTypeArray = svTypeArray;
		this.quality = quality;
	}

	public SVRecord(String chrom, int start, int end, SVType[] svTypeArray, String evidence) {
		super();
		this.chrom = chrom;
		this.start = start;
		this.end = end;
		this.svTypeArray = svTypeArray;
		this.evidence = evidence;
	}

	public SVRecord(String chrom, int start, int end, SVType[] svTypeArray, double quality, String evidence) {
		super();
		this.chrom = chrom;
		this.start = start;
		this.end = end;
		this.svTypeArray = svTypeArray;
		this.quality = quality;
		this.evidence = evidence;
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

	public int getLength() {
		return end - start + 1;
	}

	public double getQuality() {
		return quality;
	}

	public SVType[] getSVTypeArray() {
		return svTypeArray;
	}
	
	public boolean isOverlap(SVRecord svRecord) {
		if (!this.chrom.equals(svRecord.getChrom())) {
			return false;
		} else {
			if (this.getStart() > svRecord.getEnd() || this.getEnd() < svRecord.getStart()) {
				return false;
			} else
				return true;
		}
	}

	@Override
	public int compareTo(SVRecord record) {
		if (!this.getChrom().equals(record.getChrom())) {
			return this.getChrom().compareTo(record.getChrom());
		} else {
			if (this.getStart() != record.getStart()) {
				return Integer.compare(this.getStart(), record.getStart());
			} else {
				return Integer.compare(this.getEnd(), record.getEnd());
			}
		}
	}
	
	public void setEvidence(String evidence) {
		this.evidence = evidence;
	}
	
    public String getEvidence() {
		return evidence;
	}

	@Override
	public String toString() {
		if (evidence == null) {
			return chrom + "\t" + start + "\t" + end + "\t" + svTypeArray[0] + "\t" + svTypeArray[1] + "\t"
					+ svTypeArray[2] + "\t" + quality;
		} else {
			return chrom + "\t" + start + "\t" + end + "\t" + svTypeArray[0] + "\t" + svTypeArray[1] + "\t"
					+ svTypeArray[2] + "\t" + quality + "\t" + evidence;
		}
	}
}
