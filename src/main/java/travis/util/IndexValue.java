package travis.util;

public class IndexValue {
	public int index = -1;
	public double value = Double.NaN;
	public IndexValue(int i, double v) {
		index = i;
		value = v;
	}
	public double magnitude() {
		if(value < 0d)
			return -value;
		return value;
	}
}

