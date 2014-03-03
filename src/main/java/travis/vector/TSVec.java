package travis.vector;

import java.util.Arrays;
import java.util.Iterator;
import java.util.TreeMap;

import travis.util.BitHacks;
import travis.util.IndexValue;
import travis.util.TaggedIndexValue;

/**
 * SVec with "tags"
 * indices are no longer ints, but 2 ints (or a long)
 * can also be used as a matrix
 * 
 * how to handle this class...
 * i know its very useful, its common to do feature indexing with 2 (or more)
 * int keys (e.g. bigrams).
 * but, these keys are not comparable with that from a DVec, which breaks Vec
 * 
 * in general, (int,int) or long keys will never be comparable with DVec keys
 * (because you can't allocate an array that big). 2^31 gets you ~2 billion keys
 * if they're doubles, thats 16 GB of weights...
 * 
 * fundamentally, it eventually needs to get boiled down to a SVec if we are going
 * to match this up with a DVec for weights
 * 
 * TODO remove this from the Vec class/interface
 */
public class TSVec extends Vec<TSVec> {
	
	public static final int DEFAULT_CAPACITY = 16;
	
	private int[] tags;		// aka rows
	private int[] indices;	// aka cols
	private double[] values;
	private int top;
	private boolean compacted;
	
	public TSVec() { this(DEFAULT_CAPACITY); }
	
	public TSVec(int capacity) {
		tags = new int[capacity];
		indices = new int[capacity];
		values = new double[capacity];
		compacted = true;
		top = 0;
	}
	
	public TSVec(int[] tags, int[] indices, double[] values) {
		this.tags = tags;
		this.indices = indices;
		this.values = values;
		top = tags.length;
		compacted = false;
	}
	
	public int capacity() { return tags.length; }
	
	/**
	 * sort indices and consolidate duplicate entries (only for sparse vectors)
	 * @param freeExtraMem will allocate new arrays as small as possible to store tags/indices/values
	 */
	public void compact(boolean freeExtraMem) {
	
		if(compacted) return;
		
		TreeMap<Long, Double> sorted = new TreeMap<Long, Double>();
		for(int i=0; i<top; i++) {
			Long key = BitHacks.pack(tags[i], indices[i]);
			Double old = sorted.get(key);
			if(old == null) old = 0d;
			sorted.put(key, old + values[i]);
		}
		
		if(freeExtraMem) {
			int n = sorted.size();
			tags = Arrays.copyOf(tags, n);
			indices = Arrays.copyOf(indices, n);
			values = Arrays.copyOf(values, n);
		}
		
		top = 0;
		for(Long key : sorted.navigableKeySet()) {
			double val = sorted.get(key);
			if(val != 0d) {
				tags[top] = BitHacks.unpackTag(key);
				indices[top] = BitHacks.unpackIndex(key);
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
	
	private int findIndexMatching(int tag, int index) {
		compact();
		long needle = BitHacks.pack(tag, index);
		return findIndexMatching(needle, 0, tags.length);
	}
	
	private int findIndexMatching(long needle, int imin, int imax) {
		assert compacted;
		while(imin < imax) {
			int imid = (imax - imin) / 2 + imin; assert(imid < imax);
			long mid = BitHacks.pack(tags[imid], indices[imid]);
			if(mid < needle)
				imin = imid + 1;
			else
				imax = imid;
		}
		if(imax == imin) {
			long found = BitHacks.pack(tags[imin], indices[imin]);
			if(found == needle) return imin;
		}
		return -1;
	}

	@Override
	public double get(int i) { return get(0, i); }
	
	public double get(int tag, int index) {
		int idx = findIndexMatching(tag, index);
		return idx < 0 ? 0d : values[idx];
	}

	@Override
	public void add(int i, double v) { add(0, i, v); }
	
	public void add(int tag, int index, double value) {
		long ti = BitHacks.pack(tag, index);
		long tiTop = BitHacks.pack(tags[top-1], indices[top-1]);
		if(ti == tiTop) {
			values[top-1] += value;
		} else {
			if(top == capacity()) grow();
			tags[top] = tag;
			indices[top] = index;
			values[top] = value;
			compacted &= (ti > tiTop);
			top++;
		}
	}

	@Override
	public void set(int i, double v) { set(0, i, v); }
	
	public void set(int tag, int index, double value) {
		int idx = findIndexMatching(tag, index);
		if(idx < 0) add(tag, index, value);
		else values[idx] += value;
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

	@SuppressWarnings("unchecked")
	@Override
	public Iterator<IndexValue> sortedUniqNonZero() {
		compact();
		return (Iterator<IndexValue>) ((Iterator<?>) iter());
	}
	
	public Iterator<TaggedIndexValue> sortedUniqNonZeroTagged() {
		compact();
		return iter();
	}
	
	@SuppressWarnings("unchecked")
	@Override
	public Iterator<IndexValue> nonZero() {
		return (Iterator<IndexValue>) ((Iterator<?>) iter());
	}
	
	public Iterator<TaggedIndexValue> nonZeroTagged() {
		return iter();
	}
	
	private Iterator<TaggedIndexValue> iter() {
		return new Iterator<TaggedIndexValue>() {
			private int i = 0;
			private TaggedIndexValue tiv = new TaggedIndexValue(-1, -1, Double.NaN);
			@Override
			public boolean hasNext() { return i < top; }
			@Override
			public TaggedIndexValue next() {
				tiv.tag = tags[i];
				tiv.index = indices[i];
				tiv.value = values[i];
				i++;
				return tiv;
			}
			@Override
			public void remove() { throw new UnsupportedOperationException(); }
		};
	}

	@Override
	public TSVec clone() {
		int n = top;	// tags.length;
		int[] t = Arrays.copyOf(tags, n);
		int[] i = Arrays.copyOf(indices, n);
		double[] v = Arrays.copyOf(values, n);
		return new TSVec(t, i, v);
	}
}
