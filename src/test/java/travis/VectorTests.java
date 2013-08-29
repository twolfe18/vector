package travis;

import java.util.Random;
import junit.framework.Assert;
import org.junit.Test;

public class VectorTests {
	
	@Test
	public void bitHacksTests() {
		Random r = new Random(9001);
		int n = 250, k = 25;
		for(int i=0; i<n; i++) {
			
			int tag = r.nextInt(k);
			int idx = r.nextInt(k);
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
		Random r = new Random(9001);
		int t = 50;
		for(int n=1; n<=1000; n*=10) {
			for(int iter=0; iter<t; iter++) {
				for(int i=0; i<n; i++) {
					int idx = r.nextInt(size);
					double val = r.nextInt(10) - 5;
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

}
