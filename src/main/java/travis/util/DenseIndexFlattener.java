package travis.util;

import java.util.Arrays;
import java.util.Iterator;

import travis.vector.SVec;
import travis.vector.TSVec;

/*
 * it seems like there are two ways to go:
 * 1) assume a known cardinality ahead of time and use edu.jhu.hlt.parma.types.TwoPartFeatureIndexer
 * 2) scan the entire dataset, estimate cardinality, and then use something like below
 *
 */

public class DenseIndexFlattener {

	private int[] cardinalities;
	
	public DenseIndexFlattener(int numTags) {
		cardinalities = new int[numTags];
		Arrays.fill(cardinalities, 1);
	}
	
	public void observe(TSVec rawVec) {
		Iterator<TaggedIndexValue> iter = rawVec.nonZeroTagged();
		while(iter.hasNext()) {
			TaggedIndexValue tiv = iter.next();
			if(tiv.index+1 > cardinalities[tiv.tag])
				cardinalities[tiv.tag] = tiv.index+1;
		}
	}
	
	public int flatIndex(int tag, int index) {
		if(tag >= cardinalities.length)
			throw new IllegalStateException("have never observed tag=" + tag);
		int offset = 1;
		for(int i=0; i<tag; i++)
			offset *= cardinalities[tag];
		return index + offset;
	}
	
	public SVec flatten(TSVec rawVec) {
		SVec into = new SVec();
		flatten(rawVec, into);
		return into;
	}
	
	public void flatten(TSVec rawVec, SVec into) {
		Iterator<TaggedIndexValue> iter = rawVec.nonZeroTagged();
		while(iter.hasNext()) {
			TaggedIndexValue tiv = iter.next();
			int idx = flatIndex(tiv.tag, tiv.index);
			into.add(idx, tiv.value);
		}
	}
}
