package travis.vector;

import java.util.Iterator;

import travis.util.IndexValue;

/**
 * SVec with "tags"
 * indices are no longer ints, but 2 ints (or a long)
 */
public class TSVec extends Vec {
	
	private int[] tags;
	private int[] indices;
	private double[] values;
	private int top;
	private boolean compacted;

	@Override
	public double get(int i) { return get(0, i); }
	
	public double get(int tag, int index) {
		// TODO
	}

	@Override
	public void add(int i, double v) { add(0, i, v); }
	
	public void add(int tag, int index, double value) {
		// TODO
	}

	@Override
	public void set(int i, double v) { set(0, i, v); }
	
	public void set(int tag, int index, double value) {
		// TODO
	}

	@Override
	public void clear() {
		top = 0;
		compacted = true;
	}

	@Override
	public void timesEquals(double d) {
		for(int i=0; i<top; i++)
			values[i] *= d;
	}

	@Override
	public void plusEquals(double s) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public Iterator<IndexValue> nonZero() {
		// TODO Auto-generated method stub
		return null;
	}
}
