package net.imglib2.algorithm.localderivatives;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.RandomAccess;
import net.imglib2.Sampler;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.view.Views;

import org.junit.Before;
import org.junit.Test;

import Jama.Matrix;

public class HessianRandomAccessTest
{

	public static void main( final String[] args ) throws Exception
	{
		final HessianRandomAccessTest test = new HessianRandomAccessTest();
		test.setUp();

		System.out.println( "\n-------------------" );
		System.out.println( "\nCentral difference." );
		System.out.println( "\n-------------------" );

		test.testRandomAccessCD2ndOrder();
		System.out.println( "\nCentral difference, 2nd order accuracy:" );
		System.out.println( String.format( "  - error on X (~x^2): %.2e ± %.2e (n=%d).",
				test.statsErrX.mean(), test.statsErrX.std(), test.statsErrX.n() ) );
		System.out.println( String.format( "  - error on Y (~y^2): %.2e ± %.2e (n=%d).",
				test.statsErrY.mean(), test.statsErrY.std(), test.statsErrY.n() ) );
		System.out.println( String.format( "  - error on Z (Cst): %.2e ± %.2e (n=%d).",
				test.statsErrZ.mean(), test.statsErrZ.std(), test.statsErrZ.n() ) );
	}

	private Img< DoubleType > img;

	private Interval interval;

	private double a;

	private double b;

	private double c;

	private MeanStd statsErrX;

	private MeanStd statsErrY;

	private MeanStd statsErrZ;

	private long px;

	private long py;

	@Before
	public void setUp() throws Exception
	{
		img = ArrayImgs.doubles( 128, 128, 64 );

		// Center.
		px = 63;
		py = 63;

		final double theta = Math.PI / 3;
		final double ct = Math.cos( theta );
		final double st = Math.sin( theta );
		final double sx = 1.;
		final double sy = 2.;
		a = ct * ct / 2. / sx / sx + st * st / 2. / sy / sy;
		b = -Math.sin( 2. * theta ) / 4. / sx / sx + Math.sin( 2. * theta ) / 4. / sy / sy;
		c = st * st / 2. / sx / sx + ct * ct / 2. / sy / sy;

		interval = FinalInterval.createMinMax( 50, 50, 20, 70, 70, 30 );
		final Cursor< DoubleType > cursor = img.localizingCursor();
		final long[] pos = new long[ img.numDimensions() ];

		while ( cursor.hasNext() )
		{
			cursor.fwd();
			cursor.localize( pos );
			final long x = pos[ 0 ] - px;
			final long y = pos[ 1 ] - py;

			final double val = ( -a * x * x - 2 * b * x * y - c * y * y );
			cursor.get().set( val );
		}
	}

	/*
	 * Central.
	 */

	@Test
	public void testRandomAccessCD2ndOrder()
	{
		final RandomAccess< Matrix > ra = HessianRandomAccess.centralDifference( img, interval, 2 );
		test( ra, 1e-12 );
	}

	private void test( final RandomAccess< Matrix > ra, final double tolerance )
	{
		final Cursor< DoubleType > cursor = Views.interval( img, interval ).localizingCursor();

		final long[] pos = new long[ img.numDimensions() ];
		final long[] posAfter = new long[ ra.numDimensions() ];

		statsErrX = new MeanStd();
		statsErrY = new MeanStd();
		statsErrZ = new MeanStd();

		while ( cursor.hasNext() )
		{
			cursor.fwd();
			cursor.localize( pos );

			ra.setPosition( cursor );
			final Matrix matrix = ra.get();

			ra.localize( posAfter );
			assertArrayEquals( "RandomAccess position changed after access.", pos, posAfter );

			final double Hxx = matrix.get( 0, 0 );
			final double Hyy = matrix.get( 1, 1 );
			final double Hzz = matrix.get( 2, 2 );

			final double expHxx = -2. * a;
			final double expHyy = -2. * c;
			final double expHzz = 0.;

			final double dxx = Math.abs( Hxx - expHxx );
			final double dyy = Math.abs( Hyy - expHyy );
			final double dzz = Math.abs( Hzz - expHzz );

			statsErrX.add( dxx );
			statsErrY.add( dyy );
			statsErrZ.add( dzz );

			assertTrue( "Accuracy violates tolerance limit.", dxx < tolerance );
			assertTrue( "Accuracy violates tolerance limit.", dyy < tolerance );
			assertTrue( "Accuracy violates tolerance limit.", dzz < tolerance );

			final double Hxy = matrix.get( 0, 1 );
			final double Hxz = matrix.get( 0, 2 );
			final double Hyz = matrix.get( 1, 2 );

			final double expHxy = -2. * b;
			final double expHxz = 0.;
			final double expHyz = 0.;

			final double dxy = Math.abs( Hxy - expHxy );
			final double dxz = Math.abs( Hxz - expHxz );
			final double dyz = Math.abs( Hyz - expHyz );

			assertTrue( "Accuracy violates tolerance limit.", dxy < tolerance );
			assertTrue( "Accuracy violates tolerance limit.", dxz < tolerance );
			assertTrue( "Accuracy violates tolerance limit.", dyz < tolerance );
		}

		// Test copy().
		@SuppressWarnings( "rawtypes" )
		final Class< ? extends RandomAccess > class1 = ra.getClass();
		@SuppressWarnings( "rawtypes" )
		final Class< ? extends Sampler > class2 = ra.copy().getClass();
		assertSame( "The class returned by copy() is not the right one.", class1, class2 );

		// Test copyRandomAccess().
		@SuppressWarnings( "rawtypes" )
		final Class< ? extends RandomAccess > class3 = ra.copyRandomAccess().getClass();
		assertSame( "The class returned by copyRandomAccess() is not the right one.", class1, class3 );
	}

}
