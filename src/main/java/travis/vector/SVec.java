package travis.vector;

import java.util.Arrays;
import java.util.Iterator;
import java.util.TreeMap;

import travis.util.IndexValue;

public class SVec extends Vec<SVec> {
	
	public static final int DEFAULT_CAPACITY = 16;
	
	private int[] indices;
	private double[] values;
	private boolean compacted;	// are indices sorted and unique?
	private int top;			// values in `indices` and `values` with indices < top are valid
	
	public SVec() { this(DEFAULT_CAPACITY); }
	
	public SVec(int capacity) {
		indices = new int[capacity];
		values = new double[capacity];
		compacted = true;
		top = 0;
	}
	
	public SVec(int[] indices, double[] values) {
		this.indices = indices;
		this.values = values;
		top = indices.length;
		compacted = false;
	}
	
	// use this for testing the speed of getting raw arrays vs using iterator
	class UnsafeSVec extends SVec {
		public int[] getIndices() { return indices; }
		public double[] getValues() { return values; }
		public int getTop() { return top; }
	}
	
	/**
	 * sort indices and consolidate duplicate entries (only for sparse vectors)
	 * @param freeExtraMem will allocate new arrays as small as possible to store tags/indices/values
	 */
	public void compact(boolean freeExtraMem) {
	
		if(compacted) return;
		
		TreeMap<Integer, Double> sorted = new TreeMap<Integer, Double>();
		for(int i=0; i<top; i++) {
			Integer key = indices[i];
			Double old = sorted.get(key);
			if(old == null) old = 0d;
			sorted.put(key, old + values[i]);
		}
		
		if(freeExtraMem) {
			int n = sorted.size();
			indices = Arrays.copyOf(indices, n);
			values = Arrays.copyOf(values, n);
		}
		
		top = 0;
		for(Integer key : sorted.navigableKeySet()) {
			double val = sorted.get(key);
			if(val != 0d) {
				indices[top] = key;
				values[top] = val;
				top++;
			}
		}
		compacted = true;
	}
	public void compact() { compact(false); }
	
	private void grow() {
		int newSize = (int)(capacity() * 1.3d + 8d);
		indices = Arrays.copyOf(indices, newSize);
		values = Arrays.copyOf(values, newSize);
	}
	
	public int capacity() { return values.length; }
	
	private int findIndexMatching(int needle, int imin, int imax) {
		assert compacted;
		while(imin < imax) {
			int imid = (imax - imin) / 2 + imin; assert(imid < imax);
			int mid = indices[imid];
			if(mid < needle)
				imin = imid + 1;
			else
				imax = imid;
		}
		if(imax == imin) {
			int found = indices[imin];
			if(found == needle) return imin;
		}
		return -1;
	}
	
	@Override
	public double get(int i) {
		compact();
		int idx = findIndexMatching(i, 0, values.length);
		return idx < 0 ? 0d : values[idx];
	}
	
	@Override
	public void add(int i, double v) {
		if(i == indices[top-1]) {
			values[top-1] += v;
		} else {
			if(top == capacity()) grow();
			indices[top] = i;
			values[top] = v;
			compacted &= (i > indices[top-1]);
			top++;
		}
	}
	
	@Override
	public void set(int i, double v) {
		compact();
		int idx = findIndexMatching(i, 0, values.length);
		if(idx < 0) add(i, v);
		else values[idx] += v;
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
		compact();
		for(int i=0; i<top; i++)
			values[i] += s;
	}
	
	@Override
	public Iterator<IndexValue> nonZero() {
		return new Iterator<IndexValue>() {
			private int i = 0;
			private IndexValue iv = new IndexValue(-1, Double.NaN);
			@Override
			public boolean hasNext() { return i < top; }
			@Override
			public IndexValue next() {
				iv.index = indices[i];
				iv.value = values[i];
				i++;
				return iv;
			}
			@Override
			public void remove() { throw new UnsupportedOperationException(); }
		};
	}

	@Override
	public SVec clone() {
		int n = top;	// indices.length;
		int[] i = Arrays.copyOf(indices, n);
		double[] v = Arrays.copyOf(values, n);
		return new SVec(i, v);
	}

}
