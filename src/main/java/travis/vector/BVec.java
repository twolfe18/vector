package travis.vector;

import java.util.BitSet;
import java.util.Iterator;

import travis.util.IndexValue;

/**
 * bit vector backed by a java.util.BitSet
 */
public class BVec extends Vec {
	
	private BitSet values;
	
	public BVec(int dimension) {
		values = new BitSet(dimension);
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

	@Override
	public void plusEquals(double s) { throw new UnsupportedOperationException(); }

	@Override
	public Iterator<IndexValue> nonZero() {
		// TODO
	}
}
