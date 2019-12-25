package cn.edu.hit.triocnv.breakpoint;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;

/**
 *
 * @author Yongzhuang Liu
 */

public class ExtractSequence {

	private ReferenceSequenceFile referenceSequenceFile;
	private String outputFile;
	private String chrom;

	public ExtractSequence(String referenceSequenceFile, String outputFile, String chrom) throws IOException {
		this.referenceSequenceFile = new IndexedFastaSequenceFile(new File(referenceSequenceFile));
		this.outputFile = outputFile;
		this.chrom = chrom;
	}

	public void getSubSequence(String contig, long start, long end) throws IOException {
		ReferenceSequence referenceSequence = referenceSequenceFile.getSubsequenceAt(contig, start, end);
		String header = ">" + contig + ":" + start + "-" + end;
		String sequence = referenceSequence.getBaseString();
		FileWriter writer = new FileWriter(outputFile);
		writer.write(header + "\n");
		writer.write(sequence);
		writer.close();
	}
}
