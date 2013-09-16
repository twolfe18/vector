package travis.vector;

import java.util.Arrays;
import java.util.Iterator;

import travis.util.IndexValue;

public class DVec extends Vec<DVec> {
	
	private double[] values;
	
	public DVec(double[] vals) {
		values = vals;
	}
	
	public DVec(int dimension) {
		values = new double[dimension];
	}

	public int dimension() { return values.length; }
	
	@Override
	public double get(int i) { return values[i]; }

	@Override
	public void add(int i, double v) { values[i] += v; }

	@Override
	public void set(int i, double v) { values[i] = v; }

	@Override
	public void clear() { Arrays.fill(values, 0); }

	@Override
	public void timesEquals(double d) {
		for(int i=0; i<values.length; i++)
			values[i] *= d;
	}

	@Override
	public void plusEquals(double s) {
		for(int i=0; i<values.length; i++)
			values[i] += s;
	}

	@Override
	public Iterator<IndexValue> nonZero() {
		return new Iterator<IndexValue>() {
			private int i = 0;
			private IndexValue iv = new IndexValue(-1, Double.NaN);
			@Override
			public boolean hasNext() { return i < values.length; }
			@Override
			public IndexValue next() {
				iv.index = i;
				iv.value = values[i];
				i++;
				return iv;
			}
			@Override
			public void remove() { throw new UnsupportedOperationException(); }
		};
	}

	@Override
	public DVec clone() {
		return new DVec(Arrays.copyOf(values, values.length));
	}
}
