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
	
	/**
	 * product-wise: this *= (scale * other)
	 */
	public void timesEquals(Vec<?> other, double scale) {
		if(other instanceof DVec) {
			DVec d = (DVec) other;
			assert(dimension() == d.dimension());
			for(int i=0; i<values.length; i++)
				values[i] *= scale * d.values[i];
		}
		else {	// sparse
			int i = 0;
			for(IndexValue iv : other) {
				for(; i<iv.index; i++)
					values[i] = 0d;
				values[i] *= scale * iv.value;
				i++;
			}
			for(; i<values.length; i++)
				values[i] = 0;
		}
	}
	
	/**
	 * product-wise: this *= other
	 */
	public void timesEquals(Vec<?> other) { timesEquals(other, 1d); }

	@Override
	public void plusEquals(double s) {
		for(int i=0; i<values.length; i++)
			values[i] += s;
	}

	@Override
	public Iterator<IndexValue> sortedUniqNonZero() {
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
