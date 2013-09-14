package travis.vector;

import java.util.Iterator;

import travis.util.IndexValue;

public abstract class Vec {
	
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
	
	/* everything below here has an implementation based on stuff above here
	   you can always override the definitions below with more efficient versions */
	
	/**
	 * an iterator of unique nonzero index-values
	 */
	public abstract Iterator<IndexValue> nonZero();
	
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
	
	public int l0Norm() {
		int nz = 0;
		Iterator<IndexValue> iter = nonZero();
		while(iter.hasNext()) { iter.next(); nz++; }
		return nz;
	}

	public double l1Norm() {
		double s = 0d;
		Iterator<IndexValue> iter = nonZero();
		while(iter.hasNext()) {
			IndexValue iv = iter.next();
			s += iv.magnitude();
		}
		return s;
	}
	
	public double l2Norm() {
		double ss = 0d;
		Iterator<IndexValue> iter = nonZero();
		while(iter.hasNext()) {
			IndexValue iv = iter.next();
			double v = iv.value;
			ss += v * v;
		}
		return Math.sqrt(ss);
	}
	
	public double lInfNorm() {
		double mm = 0d;
		Iterator<IndexValue> iter = nonZero();
		while(iter.hasNext()) {
			IndexValue iv = iter.next();
			double m = iv.magnitude();
			if(m > mm) mm = m;
		}
		return mm;
	}
}

