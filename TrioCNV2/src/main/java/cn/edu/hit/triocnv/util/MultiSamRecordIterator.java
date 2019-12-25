package cn.edu.hit.triocnv.util;

import htsjdk.samtools.SamReader;
import java.util.Collection;
import java.util.PriorityQueue;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.util.CloseableIterator;

/**
 *
 * @author Yongzhuang Liu
 */
public class MultiSamRecordIterator {

	private final PriorityQueue<ComparableSamRecordIterator> pq;

	public MultiSamRecordIterator(Collection<SamReader> readers) {
		this.pq = new PriorityQueue<ComparableSamRecordIterator>(readers.size());
		for (SamReader reader : readers) {
			addIfNotEmpty(
					new ComparableSamRecordIterator(reader, reader.iterator(), new SAMRecordCoordinateComparator()));
		}
	}

	private void addIfNotEmpty(final ComparableSamRecordIterator iterator) {
		if (iterator.hasNext()) {
			pq.offer(iterator);
		} else {
			iterator.close();
		}
	}

	public boolean hasNext() {
		return !this.pq.isEmpty();
	}

	public SAMRecord next() {
		ComparableSamRecordIterator iterator = this.pq.poll();
		SAMRecord record = iterator.next();
		addIfNotEmpty(iterator);
		return record;
	}

	public void close() {
		for (CloseableIterator<SAMRecord> iterator : pq) {
			iterator.close();
		}
	}
}
