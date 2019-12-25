package cn.edu.hit.triocnv.util;

import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import java.util.Comparator;
import htsjdk.samtools.SAMRecord;
/**
*
* @author Yongzhuang Liu
*/
class ComparableSamRecordIterator extends PeekableIterator<SAMRecord>
		implements Comparable<ComparableSamRecordIterator> {

	private final Comparator<SAMRecord> comparator;
	private final SamReader reader;

	public ComparableSamRecordIterator(final SamReader sam, final SAMRecordIterator samRecordIterator,
			final Comparator<SAMRecord> comparator) {
		super(samRecordIterator);
		this.reader = sam;
		this.comparator = comparator;
	}

	public SamReader getReader() {
		return reader;
	}

	public int compareTo(final ComparableSamRecordIterator that) {
		if (this.comparator.getClass() != that.comparator.getClass()) {
			throw new IllegalStateException("Attempt to compare two ComparableSAMRecordIterators that "
					+ "have different orderings internally");
		}

		final SAMRecord record = this.peek();
		final SAMRecord record2 = that.peek();
		return comparator.compare(record, record2);
	}

	@Override
	public boolean equals(final Object o) {
		if (this == o) {
			return true;
		}
		if (o == null || getClass() != o.getClass()) {
			return false;
		}

		return compareTo((ComparableSamRecordIterator) o) == 0;
	}

	@Override
	public int hashCode() {
		throw new UnsupportedOperationException(
				"ComparableSamRecordIterator should not be hashed because it can change value");
	}
}
