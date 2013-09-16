package travis.vector;

import java.util.Iterator;
import java.util.List;

import travis.util.IndexValue;

/**
 * binary operators that allocate memory
 */
public class VecOps {

	/**
	 * rules (where D=DVec, S=SVec, B=BVec):
	 * D + D = D
	 * D + S = D
	 * D + B = D
	 * S + S = S
	 * S + B = S
	 * B + B = B
	 */
	public Vec<?> sum(Vec<?> a, Vec<?> b) {
		
		assert(!(a instanceof TSVec) && !(b instanceof TSVec));
		
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
		Vec<?> accum = vecs.get(0).clone();
		for(int i=1; i<vecs.size(); i++)
			accum = sum(accum, vecs.get(i));
		return accum;
	}
	
	/**
	 * rules (where D=DVec, S=SVec, B=BVec):
	 * D * D = D
	 * D * S = S
	 * D * B = S
	 * S * S = S
	 * S * B = S
	 * B * B = B
	 */
	public Vec<?> prod(Vec<?> a, Vec<?> b) {
		
		assert(!(a instanceof TSVec) && !(b instanceof TSVec));
		
		// D * D = D
		if(a instanceof DVec && b instanceof DVec)
			return ddProd((DVec) a, (DVec) b);
		
		// D * S = S
		if(a instanceof DVec && b instanceof SVec)
			return dProd((DVec) a, (SVec) b);
		if(a instanceof SVec && b instanceof DVec)
			return dProd((DVec) b, (SVec) a);
		
		// D * B = S
		if(a instanceof DVec && b instanceof BVec)
			return dProd((DVec) a, (BVec) b);
		if(a instanceof BVec && b instanceof DVec)
			return dProd((DVec) b, (BVec) a);
		
		// S * S = S
		// S * B = S
		if(a instanceof SVec)
			return sProd((SVec) a, b);
		if(b instanceof SVec)
			return sProd((SVec) b, a);
		
		// B * B = B
		return bbProd((BVec) a, (BVec) b);
	}
	
	private DVec ddProd(DVec a, DVec b) {
		DVec c = a.clone();
		c.timesEquals(b);
		return c;		
	}
	
	private SVec dProd(DVec a, Vec<?> sparse) {
		SVec c = new SVec();
		for(IndexValue iv : sparse) {
			int i = iv.index;
			double v = iv.value * a.get(i);
			c.add(i, v);
		}
		return c;
	}
	
	private SVec sProd(SVec a, Vec<?> sparse) {
		SVec result = new SVec(a.capacity());
		Iterator<IndexValue> bi = sparse.sortedUniqNonZero();
		for(IndexValue iv : a) {
			long needle = iv.fullIndex();
			IndexValue cur = null;
			while(cur != null && cur.fullIndex() < needle) {
				if(!bi.hasNext())
					return result;
				cur = bi.next();
				if(cur.fullIndex() == needle) {
					int i = iv.index;
					double v = iv.value * cur.value;
					result.add(i, v);
				}
			}
		}
		return result;
	}
	
	private BVec bbProd(BVec a, BVec b) {
		BVec c = a.clone();
		c.timesEquals(b);
		return c;
	}
	
	public Vec<?> prod(List<Vec<?>> vecs) {
		Vec<?> accum = vecs.get(0).clone();
		for(int i=1; i<vecs.size(); i++)
			accum = prod(accum, vecs.get(i));
		return accum;
	}
	
	// TODO you can also do this with a hashmap if you
	// have at least one sparse vec that is not sorted
	public double dot(Vec<?> a, Vec<?> b) {
		if(a instanceof DVec)
			return dot((DVec) a, b);
		if(b instanceof DVec)
			return dot((DVec) b, a);
		
		// general sparse impl
		double d = 0d;
		Iterator<IndexValue> bi = b.nonZero();
		for(IndexValue iv : a) {
			long needle = iv.fullIndex();
			IndexValue cur = null;
			while(cur != null && cur.fullIndex() < needle) {
				if(!bi.hasNext())
					return d;
				cur = bi.next();
				if(cur.fullIndex() == needle)
					d += iv.value * cur.value;
			}
		}
		return d;
	}
	
	// TODO benchmark this against a version that can get the backing array
	public double dot(DVec a, Vec<?> b) {
		double d = 0d;
		Iterator<IndexValue> iter = b.nonZero();
		while(iter.hasNext()) {
			IndexValue iv = iter.next();
			d += a.get(iv.index) * iv.value;
		}
		return d;
	}

}
