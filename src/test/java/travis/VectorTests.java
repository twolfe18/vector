package travis;

import java.util.Random;
import junit.framework.Assert;
import org.junit.Test;

public class VectorTests {
	
	private final Random rand = new Random(9001);
	
	@Test
	public void bitHacksTests() {
		int n = 250, k = 25;
		for(int i=0; i<n; i++) {
			
			int tag = rand.nextInt(k);
			int idx = rand.nextInt(k);
			long packed = Vector.pack(tag, idx);
			int unpackedTag = Vector.unpackTag(packed);
			int unpackedIdx = Vector.unpackIndex(packed);
			
			//System.out.printf("tag=%s idx=%s packed=%s unpackedTag=%s unpackedIdx=%s\n",
			//		Integer.toBinaryString(tag), Integer.toBinaryString(idx), Long.toBinaryString(packed),
			//		Integer.toBinaryString(unpackedTag), Integer.toBinaryString(unpackedIdx));
			
			//System.out.printf("mask=%s\n", Long.toBinaryString((1l<<32)-1));
			
			Assert.assertEquals(tag, unpackedTag);
			Assert.assertEquals(idx, unpackedIdx);
		}
	}
	
	@Test
	public void basicTests() {
		int size = 50;
		Vector d = new Vector(size);
		d.add(0, 1d);
		d.add(2, -2d);
		d.add(2, 4d);
		Assert.assertEquals(1d, d.get(0));
		Assert.assertEquals(0d, d.get(1));
		Assert.assertEquals(2d, d.get(2));
		
		d.clear();
		Vector s = new Vector(false);
		int t = 50;
		for(int n=1; n<=1000; n*=10) {
			for(int iter=0; iter<t; iter++) {
				for(int i=0; i<n; i++) {
					int idx = rand.nextInt(size);
					double val = rand.nextInt(10) - 5;
					d.add(idx, val);
					s.add(idx, val);
				}
				for(int i=0; i<size; i++)
					Assert.assertEquals(d.get(i), s.get(i));

				int dl0 = d.l0Norm();
				int sl0 = s.l0Norm();
				if(dl0 != sl0)
					System.out.printf("sparse-l0=%d dense-l0=%d\n", s.l0Norm(), d.l0Norm());

				Assert.assertEquals(d.l0Norm(), s.l0Norm());
				Assert.assertEquals(d.l2Norm(), s.l2Norm());
				Assert.assertEquals(d.l1Norm(), s.l1Norm());
				Assert.assertEquals(d.lInfNorm(), s.lInfNorm());
				
				d.clear();
				s.clear();
				Assert.assertEquals(0, d.l0Norm());
				Assert.assertEquals(0d, d.l1Norm());
				Assert.assertEquals(0d, d.l2Norm());
				Assert.assertEquals(0d, d.lInfNorm());
				Assert.assertEquals(0, s.l0Norm());
				Assert.assertEquals(0d, s.l1Norm());
				Assert.assertEquals(0d, s.l2Norm());
				Assert.assertEquals(0d, s.lInfNorm());
			}
		}
	}

	@Test
	public void normTests() {
		Vector d = new Vector(10);
		d.add(0, 1d);
		d.add(2, -2d);
		Assert.assertEquals(Math.sqrt(5d), d.l2Norm());
		Assert.assertEquals(3d, d.l1Norm());
		Assert.assertEquals(2d, d.lInfNorm());
		Assert.assertEquals(2, d.l0Norm());
		
		Vector s = new Vector(false);
		s.add(0, 1d);
		s.add(2, -2d);
		Assert.assertEquals(Math.sqrt(5d), s.l2Norm());
		Assert.assertEquals(3d, s.l1Norm());
		Assert.assertEquals(2d, s.lInfNorm());
		Assert.assertEquals(2, s.l0Norm());
		
		Vector t = new Vector(true);
		t.add(0, 0, 1d);
		t.add(1, 2, -2d);
		Assert.assertEquals(Math.sqrt(5d), t.l2Norm());
		Assert.assertEquals(3d, t.l1Norm());
		Assert.assertEquals(2d, t.lInfNorm());
		Assert.assertEquals(2, t.l0Norm());
	}
	
	@Test
	public void dotTests() {
		Vector ones = Vector.rep(1d, 10);
		Vector s = new Vector(false);
		Assert.assertEquals(0d, ones.dot(s));
		s.add(0, 1d);
		Assert.assertEquals(1d, ones.dot(s));
		s.add(2, 2d);
		Assert.assertEquals(3d, ones.dot(s));
		s.add(4, -3d);
		Assert.assertEquals(0d, ones.dot(s));
		s.add(6, 4d);
		Assert.assertEquals(4d, ones.dot(s));
		
		int tries = 50;
		for(int pow=1; pow<=3; pow++) {
			int nonzero = (int) Math.pow(10, pow), range = (int)(10d * Math.pow(10, pow));
			for(int t=0; t<tries; t++) {
				Vector d1 = new Vector(range);
				Vector d2 = new Vector(range);
				Vector s1 = new Vector(false);
				Vector s2 = new Vector(false);
				for(int i=0; i<nonzero; i++) {
					int i1 = rand.nextInt(range);
					int i2 = rand.nextInt(range);
					double v1 = rand.nextGaussian();
					double v2 = rand.nextGaussian();
					d1.add(i1, v1); s1.add(i1, v1);
					d2.add(i2, v2); s2.add(i2, v2);
				}
				double dd = d1.dot(d2);
				double sd = s1.dot(s2);
				//System.out.printf("dense=%.4f sparse=%.4f\n", dd, sd);
				Assert.assertEquals(dd, sd);
			}
		}
	}
	
	@Test
	public void probTest() {
		Vector v = Vector.rep(1d, 3);
		v.makeProbDist();
		Assert.assertEquals(1/3d, v.get(0));
		Assert.assertEquals(1/3d, v.get(1));
		Assert.assertEquals(1/3d, v.get(2));
		
		v.clear();
		//v = Vector.rep(30d, 3);
		v.makeProbDist();
		Assert.assertEquals(1/3d, v.get(0));
		Assert.assertEquals(1/3d, v.get(1));
		Assert.assertEquals(1/3d, v.get(2));
		
		v.clear();
		//v = Vector.rep(-30d, 3);
		v.makeProbDist();
		Assert.assertEquals(1/3d, v.get(0));
		Assert.assertEquals(1/3d, v.get(1));
		Assert.assertEquals(1/3d, v.get(2));
		
		v = Vector.rep(0d, 2);
		v.set(0, 1d);
		v.makeProbDist();
		double p = Math.exp(1d) / (1d + Math.exp(1d));
		assertClose(p, v.get(0));
		assertClose(1d - p, v.get(1));
		
		v.clear();
		v.add(0, 9001d);
		v.makeProbDist();
		assertClose(1d, v.get(0));
		assertClose(0d, v.get(1));
		
		v.clear();
		v.add(0, -9001d);
		v.makeProbDist();
		assertClose(0d, v.get(0));
		assertClose(1d, v.get(1));
		
		v.clear();
		v.add(0, 1d);
		v.makeProbDist(0.0001d);
		assertClose(1d, v.get(0));
		assertClose(0d, v.get(1));
		
		v.clear();
		v.add(0, 1d);
		v.makeProbDist(9999d);
		assertClose(1/2d, v.get(0), 1e-4);
		assertClose(1/2d, v.get(1), 1e-4);
	}
	
	public static void assertClose(double expected, double actual) { assertClose(expected, actual, 1e-8); }
	public static void assertClose(double expected, double actual, double epsilon) {
		double diff = actual - expected;
		Assert.assertTrue(
				String.format("expected=%.3g actual=%.3g diff=%.3g > %.3g", expected, actual, diff, epsilon),
				Math.abs(diff) <= epsilon);
	}

}
