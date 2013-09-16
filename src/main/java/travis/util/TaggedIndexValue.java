package travis.util;

public class TaggedIndexValue extends IndexValue {
	public int tag = -1;
	public TaggedIndexValue(int t, int i, double v) {
		super(i, v);
		tag = t;		
	}
	@Override
	public long fullIndex() {
		return BitHacks.pack(tag, index);
	}
}

