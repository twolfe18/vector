package travis.vector;

import java.util.BitSet;
import java.util.Iterator;

import travis.util.IndexValue;

/**
 * bit vector backed by a java.util.BitSet
 * (should have amortized cost of 1 bit per dimension)
 * 
 * TODO: should `add` have `or` semantics (for only BVec)?
 *       preliminary (conservative) answer: no
 */
public class BVec extends Vec<BVec> {
	
	private BitSet values;
	
	public BVec(int dimension) {
		values = new BitSet(dimension);
	}
	
	public BVec(BitSet values) {
		this.values = values;
	}
	
	public int dimension() {
		return values.cardinality();
	}

	@Override
	public double get(int i) {
		return values.get(i) ? 1d : 0d;
	}

	@Override
	public void add(int i, double v) {
		if(v != 1d)
			throw new IllegalArgumentException();
		assert(get(i) == 0d);
		values.set(i);
	}

	@Override
	public void set(int i, double v) {
		if(v != 1d && v != 0d)
			throw new IllegalArgumentException();
		values.set(i);
	}

	@Override
	public void clear() { values.clear(); }

	@Override
	public void timesEquals(double d) { throw new UnsupportedOperationException(); }

	public void timesEquals(BVec other) {
		values.and(other.values);
	}
	
	@Override
	public void plusEquals(double s) { throw new UnsupportedOperationException(); }

	@Override
	public Iterator<IndexValue> sortedUniqNonZero() {
		return new Iterator<IndexValue>() {
			private int shownUpTo = 0;
			private IndexValue iv = new IndexValue(-1, Double.NaN);
			@Override
			public boolean hasNext() { return shownUpTo >= 0; }
			@Override
			public IndexValue next() {
				iv.index = shownUpTo;
				iv.value = 1d;
				shownUpTo = values.nextSetBit(shownUpTo);
				return iv;
			}
			@Override
			public void remove() { throw new UnsupportedOperationException(); }
		};
	}

	@Override
	public BVec clone() {
		BitSet nbs = new BitSet(values.cardinality());
		nbs.or(values);
		return new BVec(nbs);
	}
}
