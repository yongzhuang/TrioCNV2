package cn.edu.hit.triocnv.readdepth;

/**
 *
 * @author Yongzhuang Liu
 */

public class EmissionException extends Exception {
	public EmissionException() {
		super("Emission probability is out of bounds!");
	}
}
