package travis;

import edu.jhu.prim.util.Lambda.FnIntDoubleToDouble;
import edu.jhu.prim.vector.IntDoubleVector;

/**
 * an adapter to get my Vector class to play nice with Matt's code.
 * @author travis
 */
public final class IntDoubleUnsortedVector extends Vector implements IntDoubleVector {
	
	private static final long serialVersionUID = 3455462422261002449L;
	
	public IntDoubleUnsortedVector() {
		super(false, Vector.defaultSparseInitCapacity);
	}
	
	@Override
	public double dot(double[] other) {
		Vector dense = new Vector(other);
		return dense.dot(this);
	}

	@Override
	public double dot(IntDoubleVector other) {
		double dot = 0d;
		assert this.isSparse();
		for(int i=0; i<top; i++) {
			int idx = super.idx[i];
			dot += super.vals[i] * other.get(idx);
		}
		return dot;
	}

	@Override
	public void apply(FnIntDoubleToDouble function) {
		throw reasonWhyICantImplThis();
	}

	@Override
	public void add(IntDoubleVector other) {
		throw reasonWhyICantImplThis();
	}

	@Override
	public void subtract(IntDoubleVector other) {
		throw reasonWhyICantImplThis();
	}

	@Override
	public void product(IntDoubleVector other) {
		throw reasonWhyICantImplThis();
	}
	
	private RuntimeException reasonWhyICantImplThis() {
		throw new RuntimeException("I can't implement some methods involving "
				+ "IntDoubleVector because I don't know either their dimension "
				+ "or the maximum non-zero index on which to call get()");
	}

	@Override
	public IntDoubleVector copy() {
		IntDoubleUnsortedVector n = new IntDoubleUnsortedVector();
		n.vals = super.vals.clone();
		n.tags = null;	assert super.tags == null;
		n.idx = super.idx.clone();
		n.top = super.top;
		n.compacted = super.compacted;
		n.printWarnings = super.printWarnings;
		return n;
	}

}
