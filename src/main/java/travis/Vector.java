package travis;

import java.util.Arrays;
import java.util.Iterator;
import java.util.TreeMap;

import travis.util.BitHacks;

/**
 * The one true vector!
 * 
 * Maintaining different types for sparse and dense vectors is too tedious,
 * so I'm just going to write one type that handles things right.
 * 
 * Indices must be >= 0 (I reserve the right to use negative numbers for evil).
 * 
 * The sparse usage of this vector can use two ints to index a value (e.g. a feature
 * template as the "tag" with a local "index" -- together they make a composite key).
 * The dense usage does not allow for this. Any binary ops between sparse and dense
 * instances must not use "tags" (i.e. two int indexing scheme).
 * 
 * ==== Binary Ops ====
 * Op		Sparse		Dense		Mixed
 * +		$			$			$
 * x		!			$			$			// sparse product (vector scale) requires compacting
 * dot		!			$			$			// sparse dot sparse requires set for intersection
 * 
 * ==== Unary Ops ====
 * Op		Sparse		Dense
 * norms	!			$		// sparse requires compacting
 * scale	$			$
 * shift	?			? 		// requires adding an global offset (not implemented yet)
 * get		$			$
 * set		!			$		// sparse requires compacting
 * clone	$			$
 * 
 * ==== Side-effect-free Ops ====
 * Op		Sparse		Dense		Mixed
 * sum		$			$			$
 * product	!			$			$		// no implementation for sparse-sparse products yet
 * 
 * @author Travis Wolfe <twolfe18@gmail.com>
 */
public class Vector {
	
	// TODO add ability to be a binary vector (with no support for tags)
	
	// TODO
	// 1) constructors that take arrays
	// 2) split into Vector interface and SVec, DVec that implement them
	// 3) BVec can also implement Vector based on a bitset
	// 4) have a static class that implements binary operators efficiently
	//		(this will need to inspect things like int[] idx, experiment with how to do this efficiently and safely)
	
	
	private double[] vals;
	
	// sparse only
	private int[] tags;
	private int[] idx;
	private int top;	// indices less than this are valid
	private boolean compacted;
	private int capacity() {
		assert tags == null || tags.length == idx.length;
		return idx.length;
	}
	
	public boolean printWarnings = true;
	
	/**
	 * Dense constructor
	 */
	public Vector(int dimension) {
		vals = new double[dimension];
		top = -1;
	}
	
	public static Vector dense(int dimension) { return new Vector(dimension); }
	
	/**
	 * repeat a value to make a constant vector
	 */
	public static Vector rep(double value, int times) {
		Vector v = new Vector(times);
		Arrays.fill(v.vals, value);
		return v;
	}
	
	/**
	 * Sparse constructor
	 */
	public Vector(boolean withTags, int initCapacity) {
		if(withTags) tags = new int[initCapacity];
		idx = new int[initCapacity];
		vals = new double[initCapacity];
		top = 0;
		compacted = true;
	}
	
	public static final int defaultSparseInitCapacity = 16;
	public static Vector sparse() { return new Vector(false, defaultSparseInitCapacity); }
	public static Vector sparse(int initCapacity) { return new Vector(false, initCapacity); }
	public static Vector sparse(boolean withTags) { return new Vector(withTags, defaultSparseInitCapacity); }
	public static Vector sparse(boolean withTags, int initCapacity) { return new Vector(withTags, initCapacity); }
	
	public Vector clone() {
		Vector v = new Vector(0);
		v.tags = tags == null ? null : Arrays.copyOf(tags, tags.length);
		v.idx = idx == null ? null : Arrays.copyOf(idx, idx.length);
		v.vals = Arrays.copyOf(vals, vals.length);
		v.top = top;
		v.compacted = compacted;
		return v;
	}
	
	public boolean isSparse() { return idx != null; }
	public boolean isDense() { return idx == null; }
	public boolean isTagged() { return tags != null; }
	
	/**
	 * for sparse instances, assumes tag=0
	 */
	public double get(int index) {
		if(isSparse()) return get(0, index);
		else return vals[index];
	}
	
	/**
	 * for sparse instances only
	 */
	public double get(int tag, int index) {
		if(!isSparse())
			throw new IllegalStateException("you cannot use a tag with a dense vector");
		// if we need to do an O(#non-zero) operation here anyway, might as well compact
		compact();
		int i = findIndexMatching(tag, index);
		if(i < 0) return 0d;
		else return vals[i];
	}
	
	/**
	 * @return -1 if not found
	 */
	private int findIndexMatching(int tag, int index) {		assert isSparse();
		// oh noes! implement my own binary search????
		// thanks java!
		compact();
		return findIndexMatching(tag, index, 0, top-1);
	}
	
	private int findIndexMatching(int tag, int index, int imin, int imax) {
		assert compacted;
		long needle = BitHacks.pack(tag, index);
		while(imin < imax) {
			int imid = (imax - imin) / 2 + imin; assert(imid < imax);
			long mid = BitHacks.pack(tags == null ? 0 : tags[imid], idx[imid]);
			if(mid < needle)
				imin = imid + 1;
			else
				imax = imid;
		}
		if(imax == imin) {
			long found = BitHacks.pack(tags == null ? 0 : tags[imin], idx[imin]);
			if(found == needle) return imin;
		}
		return -1;
	}
	
	private long packedIndex(int position) {
		assert isSparse();
		int tag = tags == null ? 0 : tags[position];
		return BitHacks.pack(tag, idx[position]);
	}
	
	/**
	 * sort indices and consolidate duplicate entries (only for sparse vectors)
	 * @param freeExtraMem will allocate new arrays as small as possible to store tags/indices/values
	 */
	private void compact(boolean freeExtraMem) {			assert isSparse();
	
		if(compacted) return;
		
		// TODO special case for no-tags && small-biggest-key => use array instead of treemap?
		TreeMap<Long, Double> sorted = new TreeMap<Long, Double>();
		for(int i=0; i<top; i++) {
			int tag = tags == null ? 0 : tags[i];
			long key = BitHacks.pack(tag, idx[i]);
			Double old = sorted.get(key);
			if(old == null) old = 0d;
			sorted.put(key, old + vals[i]);
		}
		
		if(freeExtraMem) {
			int n = sorted.size();
			tags = tags == null ? null : Arrays.copyOf(tags, n);
			idx = Arrays.copyOf(idx, n);
			vals = Arrays.copyOf(vals, n);
		}
		
		top = 0;
		for(Long key : sorted.navigableKeySet()) {
			double val = sorted.get(key);
			if(val != 0d) {
				if(tags != null)
					tags[top] = BitHacks.unpackTag(key);
				idx[top] = BitHacks.unpackIndex(key);
				vals[top] = val;
				top++;
			}
		}
		compacted = true;
	}
	
	private void compact() { compact(false); }
	
	/**
	 * sets this vector to the 0 vector
	 */
	public void clear() {
		if(isDense())
			Arrays.fill(vals, 0d);
		else top = 0;
		compacted = true;
	}
	
	/**
	 * for sparse instances assumes tag=0
	 */
	public double set(int index, double value) {
		if(isSparse())
			return set(0, index, value);
		else {
			double old = vals[index];
			vals[index] = value;
			return old;
		}
	}
	
	/**
	 * for sparse instances only
	 */
	public double set(int tag, int index, double value) {
		if(!isSparse())
			throw new IllegalStateException("you cannot use a tag with a dense vector");
		compact();
		int i = findIndexMatching(tag, index);
		if(i < 0) {
			add(tag, index, value);
			compacted = false;
			return 0d;
		} else {
			double old = vals[i];
			vals[i] = value;
			return old;
		}
	}
	
	/**
	 * for sparse instances, assumes tag=0
	 */
	public void add(int index, double value) {
		if(isDense())
			vals[index] += value;
		else add(0, index, value);
	}
	
	/**
	 * for sparse instances only
	 */
	public void add(int tag, int index, double value) {
		if(!isSparse())
			throw new IllegalStateException("you cannot use a tag with a dense vector");
		if(top == capacity()) grow();
		if(tags != null) tags[top] = tag;
		idx[top] = index;
		vals[top] = value;
		top++;
		compacted = false;
	}
	
	private void grow() {
		int newSize = (int)(capacity() * 1.3d + 8d);
		if(tags != null)
			tags = Arrays.copyOf(tags, newSize);
		idx = Arrays.copyOf(idx, newSize);
		vals = Arrays.copyOf(vals, newSize);
	}
	
	public int l0Norm() {
		if(isSparse()) {
			compact();
			return top;
		} else {
			int nonzero = 0;
			for(int i=0; i<vals.length; i++)
				if(vals[i] != 0d) nonzero++;
			return nonzero;
		}
	}
	
	public double l1Norm() {
		double sum = 0d;
		int n = vals.length;
		if(isSparse()) {
			compact();
			n = top;
		}
		for(int i=0; i<n; i++) {
			double v = vals[i];
			if(v >= 0d) sum += v;
			else sum -= v;
		}
		return sum;
	}
	
	public double l2Norm() {
		double sum = 0d;
		int n = vals.length;
		if(isSparse()) {
			compact();
			n = top;
		}
		for(int i=0; i<n; i++) {
			double v = vals[i];
			sum += v * v;
		}
		return Math.sqrt(sum);
	}
	
	public double lInfNorm() {
		double biggest = 0d;
		int n = vals.length;
		if(isSparse()) {
			compact();
			n = top;
		}
		for(int i=0; i<n; i++) {
			double v = vals[i];
			if(v < 0d && v < biggest)
				biggest = v;
			else if(v > 0d && v > biggest)
				biggest = v;
		}
		return biggest >= 0 ? biggest : -biggest;
	}
	
	public void makeUnitVector() { scale(1d / l2Norm()); }
	
	/**
	 * exponentiate and re-normalize (only for dense vectors)
	 * temperature = 0 is Kronecker delta on highest value
	 * temperature = 1 is regular exp & re-normalize
	 * temperature = infinity is uniform distribution
	 */
	public void makeProbDist(double temperature) {
		if(temperature < 1e-8)
			throw new IllegalArgumentException("that temperature is too low!");
		if(isSparse())
			throw new RuntimeException("dimension sparse vec not specified, can't compute Z");
		
		// exp(x) = exp(B) * exp(x - B)
		// apply temperature change first, then apply this^ rule
		double max = vals[0];
		for(int i=1; i<vals.length; i++) {
			double v = vals[i];
			if(v > max) max = v;
		}
		double B = max / temperature;
		
		double Z = 0d;
		for(int i=0; i<vals.length; i++) {
			double v = Math.exp(vals[i] / temperature - B);
			vals[i] = v;
			Z += v;
		}
		scale(1d / Z);
	}
	
	public void makeProbDist() { makeProbDist(1d); }
	
	public void exponentiate(double power) {
		int n = vals.length;
		if(isSparse()) {
			n = top;
			compact();
		}
		for(int i=0; i<n; i++) vals[i] = Math.pow(vals[i], power);
	}
	
	public void scale(double factor) {
		// no need to compact here: a*x + a*y = a*(x+y)
		int n = isSparse() ? top : vals.length;
		for(int i=0; i<n; i++) vals[i] *= factor;
	}
	
	/**
	 * component-wise multiplication
	 */
	public void scale(Vector other) {
		if(isSparse() && other.isSparse()) {
			// use set intersection
			// (maybe compact both vecs?)
			throw new RuntimeException("implement me");
		}
		else if(isSparse() && other.isDense()) {
			if(tags != null)
				throw new RuntimeException("can't coerce sparse vec with tags to dense vec");
			// no need to compact here: a*x + a*y = a*(x+y)
			for(int i=0; i<top; i++)
				vals[i] *= other.vals[idx[i]];
		}
		else if(isDense() && other.isSparse()) {
			// no need to compact here: a*x + a*y = a*(x+y)
			if(other.tags != null)
				throw new RuntimeException("can't coerce sparse vec with tags to dense vec");
			int ptr = 0;
			for(int i=0; i<other.top; i++) {
				int idx = other.idx[i];
				while(ptr < idx) vals[ptr++] = 0d;
				vals[ptr++] *= other.vals[i];
			}
			while(ptr < vals.length) vals[ptr++] = 0d;
		}
		else {
			if(vals.length != other.vals.length)
				throw new RuntimeException("dimensions must match");
			for(int i=0; i<vals.length; i++)
				vals[i] *= other.vals[i];
		}
	}
	
	
	/**
	 * this += scale * other
	 */
	public void add(Vector other, double scale) {
		if(isSparse() && other.isSparse()) {
			for(int i=0; i<other.top; i++) {
				int t = other.tags == null ? 0 : other.tags[i];
				add(t, other.idx[i], scale * other.vals[i]);
			}
		}
		else if(isSparse() && other.isDense()) {
			if(printWarnings)
				System.err.println("WARNING: `sparse += dense` is highly inefficient");
			if(tags != null)
				throw new RuntimeException("can't coerce sparse vec with tags to dense vec");
			for(int i=0; i<other.vals.length; i++)
				add(i, other.get(i));
		}
		else if(isDense() && other.isSparse()) {
			if(other.tags != null)
				throw new RuntimeException("can't coerce sparse vec with tags to dense vec");
			for(int i=0; i<other.top; i++)
				add(other.idx[i], scale * other.vals[i]);
		}
		else {
			if(vals.length != other.vals.length)
				throw new RuntimeException("dimensions must match");
			for(int i=0; i<vals.length; i++)
				vals[i] += scale * other.vals[i];
		}
	}
	
	public void add(Vector other) { add(other, 1d); }
	
	public double dot(Vector other) {
		double dot = 0d;
		if(isSparse() && other.isSparse()) {
			Vector smaller = this, bigger = other;
			if(this.top > other.top) {
				smaller = other; bigger = this;
			}
			smaller.compact();
			bigger.compact();
			//System.out.printf("smaller.top=%d bigger.top=%d\n", smaller.top, bigger.top);
			int j = 0;
			long attempt = bigger.packedIndex(j);
			for(int i=0; i<smaller.top; i++) {
				long needle = smaller.packedIndex(i);
				while(attempt < needle && j < bigger.top-1)
					attempt = bigger.packedIndex(++j);
				if(attempt == needle)
					dot += smaller.vals[i] * bigger.vals[j];
				if(j == bigger.top)
					break;
			}
		}
		else if(isSparse() && other.isDense()) {
			dot = other.dot(this);
		}
		else if(isDense() && other.isSparse()) {
			// no need to compact here: a*x + a*y = a*(x+y)
			if(other.tags != null)
				throw new RuntimeException("can't coerce sparse vec with tags to dense vec");
			for(int i=0; i<other.top; i++)
				dot += get(other.idx[i]) * other.vals[i];
		}
		else {
			if(vals.length != other.vals.length)
				throw new RuntimeException("dimensions must match");
			for(int i=0; i<vals.length; i++)
				dot += vals[i] * other.vals[i];
		}
		return dot;
	}
	
	/**
	 * component-wise product
	 */
	public static Vector prod(Vector a, Vector b) {
		if(a.isSparse() /*&& b.isDense()*/); { Vector t = a; a = b; b = t; }
		Vector c = a.clone();
		c.scale(b);
		return c;
	}

	/**
	 * scalar-vector multiplication
	 */
	public static Vector prod(Vector a, double factor) {
		Vector b = a.clone();
		b.scale(factor);
		return b;
	}

	/**
	 * allocates a new vector, doesn't modify arguments
	 */
	public static Vector sum(Vector a, double aScale, Vector b, double bScale) {
		if(a.isSparse()) {
			Vector t = a; a = b; b = t;
			double tt = aScale; aScale = bScale; bScale = tt;
		}
		Vector c = a.clone();
		c.scale(aScale);
		c.add(b, bScale);
		return c;
	}
	
	public static Vector sum(Vector a, Vector b) {
		return sum(a, 1d, b, 1d);
	}
	
	/**
	 * a and b must not have tags (but can be sparse or dense)
	 */
	public static Vector outerProduct(Vector a, Vector b) {
		int na = a.isSparse() ? a.top : a.vals.length;
		int nb = b.isSparse() ? b.top : b.vals.length;
		Vector s = Vector.sparse(true, na * nb);
		outerProduct(a, b, s);
		return s;
	}
	
	/**
	 * result must be a sparse and tagged.
	 * adds to result, which must be 0 for correctness.
	 * indices from a will go in the tag, indices from b will go in
	 * the index (of result).
	 */
	public static void outerProduct(Vector a, Vector b, Vector result) {
		if(!(result.isSparse() && result.isTagged()))
			throw new IllegalArgumentException("must provide a sparse tagged result vector");
		boolean aSparse = a.isSparse(), bSparse = b.isSparse();
		int na = aSparse ? a.top : a.vals.length;
		int nb = bSparse ? b.top : b.vals.length;
		for(int i=0; i<na; i++) {
			int idx_i = aSparse ? a.idx[i] : i;
			double val_i = a.vals[i];
			for(int j=0; j<nb; j++) {
				int idx_j = bSparse ? b.idx[j] : j;
				double value = val_i * b.vals[j];
				result.add(idx_i, idx_j, value);
			}
		}
	}
	
	public Vector idxPairs() {
		int n = isSparse() ? top : vals.length;
		Vector pairs = Vector.sparse(false, n * (n-1) / 2);
		idxPairs(pairs);
		return pairs;
	}
	
	public void idxPairs(Vector addTo) {
		boolean s = isSparse();
		int n = s ? top : vals.length;
		for(int i=0; i<n-1; i++) {
			int ii = s ? idx[i] : i;
			double iv = vals[i];
			for(int j=i+1; j<n; j++) {
				int jj = s ? idx[j] : j;
				double jv = vals[j];
				// index int overflow is OK because we
				// are going to hash these values anyway
				addTo.add(ii * jj, iv * jv);
			}
		}
	}
	
	private class DenseIdxIter implements Iterator<Integer> {
		private int dim, cur = 0;
		public DenseIdxIter(int dimension) { dim = dimension; }
		@Override
		public boolean hasNext() { return cur < dim; }
		@Override
		public Integer next() { return cur++; }
		@Override
		public void remove() { throw new UnsupportedOperationException(); }
	}
	
	private class SparseIdxIter implements Iterator<Integer> {
		private int i = 0, top;
		private int[] idx;
		public SparseIdxIter(int[] idx, int top) {
			this.idx = idx;
			this.top = top;
		}
		@Override
		public boolean hasNext() { return i < top; }
		@Override
		public Integer next() { return idx[i++]; }
		@Override
		public void remove() { throw new UnsupportedOperationException(); }
	}
	
	public Iterator<Integer> indices() {
		if(isTagged()) throw new RuntimeException();
		if(isSparse()) return new SparseIdxIter(idx, top);
		else return new DenseIdxIter(vals.length);
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		if(isSparse()) {
			sb.append("{");
			for(int i=0; i<top; i++) {
				String tag = tags == null ? "" : tags[i] + ":";
				sb.append(String.format("%s%d:%.2f", tag, idx[i], vals[i]));
				if(i < top - 1) sb.append(", ");
			}
			sb.append("}");
		}
		else {
			sb.append("[");
			for(int i=0; i<vals.length; i++) {
				sb.append(String.format("%.2f", vals[i]));
				if(i < vals.length - 1) sb.append(", ");
			}
			sb.append("]");
		}
		return sb.toString();
	}

	/**
	 * returns true if any values are NaN or Inf
	 */
	public boolean hasBadValues() {
		int n = isSparse() ? top : vals.length;
		for(int i=0; i<n; i++) {
			double v = vals[i];
			boolean bad = Double.isNaN(v) || Double.isInfinite(v);
			if(bad) return true;
		}
		return false;
	}

}
