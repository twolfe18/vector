package travis.util;

public class BitHacks {

	public static long pack(int tag, int index) { return ((long)tag << 32) | index; }
	public static int unpackTag(long key) { return (int)(key >>> 32); }
	public static int unpackIndex(long key) { return (int)(key & ((1l<<32)-1l)); }
	
}
