package travis.vector;

import java.util.Iterator;

import travis.util.IndexValue;

/**
 * @param <T> will be either DVec, SVec, BVec, or TSVec
 */
public abstract class Vec<T extends Vec<?>> implements Iterable<IndexValue> {
	
	public abstract double get(int i);
	
	public abstract void add(int i, double v);
	
	/**
	 * use add instead if at all possible
	 * (which is more efficient for sparse)
	 */
	public abstract void set(int i, double v);
	
	/**
	 * set all values to 0
	 */
	public abstract void clear();
	
	public abstract void timesEquals(double d);
	
	public abstract void plusEquals(double s);
	
	/**
	 * deep copy
	 */
	public abstract T clone();
	
	/**
	 * an iterator of unique nonzero index-values where the
	 * fullIndex()s returned are _unique_ and _sorted_ (ascending)
	 * 
	 * NOTE: you should capture (and use) references to
	 * things returned from calls to next() (these iterators
	 * may avoid allocating a new IndexValue for performance).
	 */
	public abstract Iterator<IndexValue> sortedUniqNonZero();

	// NOTE I could introduce:
	// public abstract Iterator<IndexValue> uniqNonZero();
	// (which is all that is needed to l0, l1, l2 norms), but
	// [1] I don't see an efficient impl of this (that is not also sortedUniqNonZero)
	// [2] if this becomes an issue, someone can override the norm functions
	
	/* everything below here has an implementation based on stuff above here
	   you can always override the definitions below with more efficient versions */
	
	/**
	 * an iterator of unique nonzero index-values
	 * 
	 * NOTE: you should capture (and use) references to
	 * things returned from calls to next() (these iterators
	 * may avoid allocating a new IndexValue for performance).
	 */
	public Iterator<IndexValue> nonZero() {
		return sortedUniqNonZero();
	}
	
	/**
	 * for use in foreach loops, calls sortedUniqNonZero()
	 */
	public final Iterator<IndexValue> iterator() {
		return sortedUniqNonZero();
	}
	
	public void makeUnit() { timesEquals(1d / l2Norm()); }
	
	public boolean hasBadValues() {
		Iterator<IndexValue> iter = nonZero();
		while(iter.hasNext()) {
			IndexValue iv = iter.next();
			if(Double.isInfinite(iv.value) || Double.isNaN(iv.value))
				return true;
		}
		return false;
	}
	
	/**
	 * this += scale * other
	 */
	public void add(Vec<?> other, double scale) {
		Iterator<IndexValue> iter = nonZero();
		while(iter.hasNext()) {
			IndexValue iv = iter.next();
			add(iv.index, scale * iv.value);
		}
	}
	
	/**
	 * this += other
	 */
	public void add(Vec<?> other) { add(other, 1d); }
	
	public int l0Norm() {
		int nz = 0;
		Iterator<IndexValue> iter = sortedUniqNonZero();
		while(iter.hasNext()) { iter.next(); nz++; }
		return nz;
	}

	public double l1Norm() {
		double s = 0d;
		Iterator<IndexValue> iter = sortedUniqNonZero();
		while(iter.hasNext()) {
			IndexValue iv = iter.next();
			s += iv.magnitude();
		}
		return s;
	}
	
	public double l2Norm() {
		double ss = 0d;
		Iterator<IndexValue> iter = sortedUniqNonZero();
		while(iter.hasNext()) {
			IndexValue iv = iter.next();
			double v = iv.value;
			ss += v * v;
		}
		return Math.sqrt(ss);
	}
	
	public double lInfNorm() {
		double mm = 0d;
		Iterator<IndexValue> iter = sortedUniqNonZero();
		while(iter.hasNext()) {
			IndexValue iv = iter.next();
			double m = iv.magnitude();
			if(m > mm) mm = m;
		}
		return mm;
	}
}

