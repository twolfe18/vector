package travis.vector;

import java.util.Iterator;
import java.util.List;

import travis.util.IndexValue;

/**
 * binary operators that allocate memory
 */
public class VecOps {

	public Vec<?> sum(Vec a, Vec b) {
		
		// D + D = D
		// D + S = D
		// D + B = D
		if(a instanceof DVec)
			return sum((DVec) a, b);
		if(b instanceof DVec)
			return sum((DVec) b, a);
		
		// S + S = S
		// S + B = S
		if(a instanceof SVec)
			return sum((SVec) a, b);
		if(b instanceof SVec)
			return sum((SVec) b, a);
		
		// B + B = B
		Vec<?> c = a.clone();
		c.add(b);
		return c;
	}
	
	/**
	 * @param a
	 * @param b: either SVec or BVec
	 */
	private DVec sum(DVec a, Vec<?> b) {
		DVec c = a.clone();
		c.add(b);
		return c;
	}
	
	/**
	 * @param a
	 * @param b: either SVec or BVec
	 * @return
	 */
	private SVec sum(SVec a, Vec<?> b) {
		SVec c = a.clone();
		c.add(b);
		return c;
	}
	
	public Vec<?> sum(List<Vec<?>> vecs) {
		Vec<?> accum = vecs.get(0);
		for(int i=1; i<vecs.size(); i++)
			accum = sum(accum, vecs.get(i));
		return accum;
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

}
