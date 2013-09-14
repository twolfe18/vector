package travis.vector;

import java.util.Iterator;

import travis.util.IndexValue;

/**
 * binary operators that allocate memory
 */
public class VecOps {

	public Vec sum(Vec a, Vec b) {
		
	}
	
	public Vec sum(List<Vec> vecs) {
		
	}
	
	public Vec prod(Vec a, Vec b) {
		
	}
	
	public Vec prod(List<Vec> vecs) {
		
	}
	
	public double dot(Vec a, Vec b) {
		
	}
	
	// TODO benchmark this against a version that can get the backing array
	public double dot(DVec a, SVec b) {
		double d = 0d;
		Iterator<IndexValue> iter = b.nonZero();
		while(iter.hasNext()) {
			IndexValue iv = iter.next();
			d += a.get(iv.index) * iv.value;
		}
		return d;
	}

	/*
	 * would be nice to implement all binary ops as:
	 * 
	 * def op(a, b):
	 *   c = a.clone
	 *   c op= b
	 *   return c
	 * 
	 * this does not force you to support op= in all derived classes!
	 */	
}
