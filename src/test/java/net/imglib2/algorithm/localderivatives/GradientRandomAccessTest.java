package net.imglib2.algorithm.localderivatives;

import static org.junit.Assert.assertTrue;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.RandomAccess;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.view.Views;

import org.junit.Before;
import org.junit.Test;

import Jama.Matrix;

public class GradientRandomAccessTest
{

	public static void main( final String[] args ) throws Exception
	{
		final GradientRandomAccessTest test = new GradientRandomAccessTest();

		test.setUp();
		test.testRandomAccess2ndOrder();
		System.out.println( "2nd order accuracy:" );
		System.out.println( String.format( "  - error on X (~1/x): %.2e ± %.2e (n=%d).",
				test.statsErrX.mean(), test.statsErrX.std(), test.statsErrX.n() ) );
		System.out.println( String.format( "  - error on Y (Cst): %.2e ± %.2e (n=%d).",
				test.statsErrY.mean(), test.statsErrY.std(), test.statsErrY.n() ) );
		System.out.println( String.format( "  - error on Z (~z): %.2e ± %.2e (n=%d).",
				test.statsErrZ.mean(), test.statsErrZ.std(), test.statsErrZ.n() ) );

		test.testRandomAccess4thOrder();
		System.out.println( "4th order accuracy:" );
		System.out.println( String.format( "  - error on X (~1/x): %.2e ± %.2e (n=%d).",
				test.statsErrX.mean(), test.statsErrX.std(), test.statsErrX.n() ) );
		System.out.println( String.format( "  - error on Y (Cst): %.2e ± %.2e (n=%d).",
				test.statsErrY.mean(), test.statsErrY.std(), test.statsErrY.n() ) );
		System.out.println( String.format( "  - error on Z (~z): %.2e ± %.2e (n=%d).",
				test.statsErrZ.mean(), test.statsErrZ.std(), test.statsErrZ.n() ) );

		test.testRandomAccess6thOrder();
		System.out.println( "6th order accuracy:" );
		System.out.println( String.format( "  - error on X (~1/x): %.2e ± %.2e (n=%d).",
				test.statsErrX.mean(), test.statsErrX.std(), test.statsErrX.n() ) );
		System.out.println( String.format( "  - error on Y (Cst): %.2e ± %.2e (n=%d).",
				test.statsErrY.mean(), test.statsErrY.std(), test.statsErrY.n() ) );
		System.out.println( String.format( "  - error on Z (~z): %.2e ± %.2e (n=%d).",
				test.statsErrZ.mean(), test.statsErrZ.std(), test.statsErrZ.n() ) );

		test.testRandomAccess8thOrder();
		System.out.println( "8th order accuracy:" );
		System.out.println( String.format( "  - error on X (~1/x): %.2e ± %.2e (n=%d).",
				test.statsErrX.mean(), test.statsErrX.std(), test.statsErrX.n() ) );
		System.out.println( String.format( "  - error on Y (Cst): %.2e ± %.2e (n=%d).",
				test.statsErrY.mean(), test.statsErrY.std(), test.statsErrY.n() ) );
		System.out.println( String.format( "  - error on Z (~z): %.2e ± %.2e (n=%d).",
				test.statsErrZ.mean(), test.statsErrZ.std(), test.statsErrZ.n() ) );
	}

	private Img< DoubleType > img;

	private F fx;

	private F fy;

	private F fz;

	private Interval interval;

	private MeanStd statsErrX;

	private MeanStd statsErrY;

	private MeanStd statsErrZ;

	@Before
	public void setUp() throws Exception
	{
		img = ArrayImgs.doubles( 128, 128, 64 );
		interval = FinalInterval.createMinMax( 50, 50, 20, 70, 70, 30 );
		final Cursor< DoubleType > cursor = img.localizingCursor();
		final long[] pos = new long[ img.numDimensions() ];

		// in X: 1/X
		final double A = 100.;
		fx = new F()
		{
			@Override
			public double eval( final double in )
			{
				return 1 / ( in + 1. );
			}
		};
		fy = new F()
		{
			@Override
			public double eval( final double in )
			{
				return A;
			}
		};
		fz = new F()
		{
			@Override
			public double eval( final double in )
			{
				return in;
			}
		};

		while ( cursor.hasNext() )
		{
			cursor.fwd();
			cursor.localize( pos );
			final long x = pos[ 0 ];
			final long y = pos[ 1 ];
			final long z = pos[ 2 ];

			final double val = fx.eval( x ) * fy.eval( y ) * fz.eval( z );
			cursor.get().set( val );
		}
	}

	@Test
	public void testRandomAccess2ndOrder()
	{
		final RandomAccess< Matrix > ra = GradientRandomAccess.centralDifference( img, interval, 2 );
		// Putting a limit on tolerance only makes sense for slowly varying function (from one pixel to the next).
		test( ra, 1e-3 );
	}

	@Test
	public void testRandomAccess4thOrder()
	{
		final RandomAccess< Matrix > ra = GradientRandomAccess.centralDifference( img, interval, 4 );
		test( ra, 1e-6 );
	}

	@Test
	public void testRandomAccess6thOrder()
	{
		final RandomAccess< Matrix > ra = GradientRandomAccess.centralDifference( img, interval, 6 );
		test( ra, 1e-8 );
	}

	@Test
	public void testRandomAccess8thOrder()
	{
		final RandomAccess< Matrix > ra = GradientRandomAccess.centralDifference( img, interval, 8 );
		test( ra, 1e-10 );
	}

	private void test( final RandomAccess< Matrix > ra, final double tolerance )
	{
		final Cursor< DoubleType > cursor = Views.interval( img, interval ).localizingCursor();

		final long[] pos = new long[ img.numDimensions() ];

		statsErrX = new MeanStd();
		statsErrY = new MeanStd();
		statsErrZ = new MeanStd();

		while ( cursor.hasNext() )
		{
			cursor.fwd();
			cursor.localize( pos );
			final long x = pos[ 0 ];
			final long y = pos[ 1 ];
			final long z = pos[ 2 ];

			ra.setPosition( cursor );
			final Matrix matrix = ra.get();
			final double gx = matrix.get( 0, 0 );
			final double gy = matrix.get( 1, 0 );
			final double gz = matrix.get( 2, 0 );
			
			final double expgx = fy.eval( y ) * fz.eval( z ) * ( -1. / ( x + 1. ) / ( x + 1. ) );
			final double expgy = 0.;
			final double expgz = fx.eval( x ) * fy.eval( y );

			final double dx = Math.abs( gx - expgx );
			final double dy = Math.abs( gy - expgy );
			final double dz = Math.abs( gz - expgz );

			statsErrX.add( dx );
			statsErrY.add( dy );
			statsErrZ.add( dz );

			assertTrue( "Accuracy violates tolerance limit.", dx < tolerance );
			assertTrue( "Accuracy violates tolerance limit.", dy < tolerance );
			assertTrue( "Accuracy violates tolerance limit.", dz < tolerance );
		}
	}

	private interface F
	{
		public double eval( double in );
	}

	private class MeanStd
	{

		private long n = 0;

		private double mean = 0.;

		private double M2 = 0.;

		public void add( final double x )
		{
			n++;
			final double delta = x - mean;
			mean += delta / n;
			M2 += delta * ( x - mean );
		}

		public double mean()
		{
			return mean;
		}

		public double std()
		{
			if ( n < 2 ) { return Double.NaN; }
			return M2 / ( n - 1 );
		}

		public long n()
		{
			return n;
		}
	}

}
