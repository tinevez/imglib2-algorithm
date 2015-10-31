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

public class GradientRandomAccessTest
{

	public static void main( final String[] args ) throws Exception
	{
		final GradientRandomAccessTest test = new GradientRandomAccessTest();
		test.setUp();

		System.out.println( "\n-------------------" );
		System.out.println( "\nCentral difference." );
		System.out.println( "\n-------------------" );

		test.testRandomAccessCD2ndOrder();
		System.out.println( "\nCentral difference, 2nd order accuracy:" );
		System.out.println( String.format( "  - error on X (~1/x): %.2e ± %.2e (n=%d).",
				test.statsErrX.mean(), test.statsErrX.std(), test.statsErrX.n() ) );
		System.out.println( String.format( "  - error on Y (Cst): %.2e ± %.2e (n=%d).",
				test.statsErrY.mean(), test.statsErrY.std(), test.statsErrY.n() ) );
		System.out.println( String.format( "  - error on Z (~z): %.2e ± %.2e (n=%d).",
				test.statsErrZ.mean(), test.statsErrZ.std(), test.statsErrZ.n() ) );

		test.testRandomAccessCD4thOrder();
		System.out.println( "\nCentral difference, 4th order accuracy:" );
		System.out.println( String.format( "  - error on X (~1/x): %.2e ± %.2e (n=%d).",
				test.statsErrX.mean(), test.statsErrX.std(), test.statsErrX.n() ) );
		System.out.println( String.format( "  - error on Y (Cst): %.2e ± %.2e (n=%d).",
				test.statsErrY.mean(), test.statsErrY.std(), test.statsErrY.n() ) );
		System.out.println( String.format( "  - error on Z (~z): %.2e ± %.2e (n=%d).",
				test.statsErrZ.mean(), test.statsErrZ.std(), test.statsErrZ.n() ) );

		test.testRandomAccessCD6thOrder();
		System.out.println( "\nCentral difference, 6th order accuracy:" );
		System.out.println( String.format( "  - error on X (~1/x): %.2e ± %.2e (n=%d).",
				test.statsErrX.mean(), test.statsErrX.std(), test.statsErrX.n() ) );
		System.out.println( String.format( "  - error on Y (Cst): %.2e ± %.2e (n=%d).",
				test.statsErrY.mean(), test.statsErrY.std(), test.statsErrY.n() ) );
		System.out.println( String.format( "  - error on Z (~z): %.2e ± %.2e (n=%d).",
				test.statsErrZ.mean(), test.statsErrZ.std(), test.statsErrZ.n() ) );

		test.testRandomAccessCD8thOrder();
		System.out.println( "\nCentral difference, 8th order accuracy:" );
		System.out.println( String.format( "  - error on X (~1/x): %.2e ± %.2e (n=%d).",
				test.statsErrX.mean(), test.statsErrX.std(), test.statsErrX.n() ) );
		System.out.println( String.format( "  - error on Y (Cst): %.2e ± %.2e (n=%d).",
				test.statsErrY.mean(), test.statsErrY.std(), test.statsErrY.n() ) );
		System.out.println( String.format( "  - error on Z (~z): %.2e ± %.2e (n=%d).",
				test.statsErrZ.mean(), test.statsErrZ.std(), test.statsErrZ.n() ) );

		System.out.println( "\n-------------------" );
		System.out.println( "\nForward difference." );
		System.out.println( "\n-------------------" );

		test.testRandomAccessFD1stOrder();
		System.out.println( "\nForward difference, 1st order accuracy:" );
		System.out.println( String.format( "  - error on X (~1/x): %.2e ± %.2e (n=%d).",
				test.statsErrX.mean(), test.statsErrX.std(), test.statsErrX.n() ) );
		System.out.println( String.format( "  - error on Y (Cst): %.2e ± %.2e (n=%d).",
				test.statsErrY.mean(), test.statsErrY.std(), test.statsErrY.n() ) );
		System.out.println( String.format( "  - error on Z (~z): %.2e ± %.2e (n=%d).",
				test.statsErrZ.mean(), test.statsErrZ.std(), test.statsErrZ.n() ) );

		test.testRandomAccessFD2ndOrder();
		System.out.println( "\nForward difference, 2nd order accuracy:" );
		System.out.println( String.format( "  - error on X (~1/x): %.2e ± %.2e (n=%d).",
				test.statsErrX.mean(), test.statsErrX.std(), test.statsErrX.n() ) );
		System.out.println( String.format( "  - error on Y (Cst): %.2e ± %.2e (n=%d).",
				test.statsErrY.mean(), test.statsErrY.std(), test.statsErrY.n() ) );
		System.out.println( String.format( "  - error on Z (~z): %.2e ± %.2e (n=%d).",
				test.statsErrZ.mean(), test.statsErrZ.std(), test.statsErrZ.n() ) );

		test.testRandomAccessFD3rdOrder();
		System.out.println( "\nForward difference, 3rd order accuracy:" );
		System.out.println( String.format( "  - error on X (~1/x): %.2e ± %.2e (n=%d).",
				test.statsErrX.mean(), test.statsErrX.std(), test.statsErrX.n() ) );
		System.out.println( String.format( "  - error on Y (Cst): %.2e ± %.2e (n=%d).",
				test.statsErrY.mean(), test.statsErrY.std(), test.statsErrY.n() ) );
		System.out.println( String.format( "  - error on Z (~z): %.2e ± %.2e (n=%d).",
				test.statsErrZ.mean(), test.statsErrZ.std(), test.statsErrZ.n() ) );

		test.testRandomAccessFD4thOrder();
		System.out.println( "\nForward difference, 4th order accuracy:" );
		System.out.println( String.format( "  - error on X (~1/x): %.2e ± %.2e (n=%d).",
				test.statsErrX.mean(), test.statsErrX.std(), test.statsErrX.n() ) );
		System.out.println( String.format( "  - error on Y (Cst): %.2e ± %.2e (n=%d).",
				test.statsErrY.mean(), test.statsErrY.std(), test.statsErrY.n() ) );
		System.out.println( String.format( "  - error on Z (~z): %.2e ± %.2e (n=%d).",
				test.statsErrZ.mean(), test.statsErrZ.std(), test.statsErrZ.n() ) );

		test.testRandomAccessFD5thOrder();
		System.out.println( "\nForward difference, 5th order accuracy:" );
		System.out.println( String.format( "  - error on X (~1/x): %.2e ± %.2e (n=%d).",
				test.statsErrX.mean(), test.statsErrX.std(), test.statsErrX.n() ) );
		System.out.println( String.format( "  - error on Y (Cst): %.2e ± %.2e (n=%d).",
				test.statsErrY.mean(), test.statsErrY.std(), test.statsErrY.n() ) );
		System.out.println( String.format( "  - error on Z (~z): %.2e ± %.2e (n=%d).",
				test.statsErrZ.mean(), test.statsErrZ.std(), test.statsErrZ.n() ) );

		test.testRandomAccessFD6thOrder();
		System.out.println( "\nForward difference, 6th order accuracy:" );
		System.out.println( String.format( "  - error on X (~1/x): %.2e ± %.2e (n=%d).",
				test.statsErrX.mean(), test.statsErrX.std(), test.statsErrX.n() ) );
		System.out.println( String.format( "  - error on Y (Cst): %.2e ± %.2e (n=%d).",
				test.statsErrY.mean(), test.statsErrY.std(), test.statsErrY.n() ) );
		System.out.println( String.format( "  - error on Z (~z): %.2e ± %.2e (n=%d).",
				test.statsErrZ.mean(), test.statsErrZ.std(), test.statsErrZ.n() ) );

		System.out.println( "\n--------------------" );
		System.out.println( "\nBackward difference." );
		System.out.println( "\n--------------------" );

		test.testRandomAccessBD1stOrder();
		System.out.println( "\nBackward difference, 1st order accuracy:" );
		System.out.println( String.format( "  - error on X (~1/x): %.2e ± %.2e (n=%d).",
				test.statsErrX.mean(), test.statsErrX.std(), test.statsErrX.n() ) );
		System.out.println( String.format( "  - error on Y (Cst): %.2e ± %.2e (n=%d).",
				test.statsErrY.mean(), test.statsErrY.std(), test.statsErrY.n() ) );
		System.out.println( String.format( "  - error on Z (~z): %.2e ± %.2e (n=%d).",
				test.statsErrZ.mean(), test.statsErrZ.std(), test.statsErrZ.n() ) );

		test.testRandomAccessBD2ndOrder();
		System.out.println( "\nBackward difference, 2nd order accuracy:" );
		System.out.println( String.format( "  - error on X (~1/x): %.2e ± %.2e (n=%d).",
				test.statsErrX.mean(), test.statsErrX.std(), test.statsErrX.n() ) );
		System.out.println( String.format( "  - error on Y (Cst): %.2e ± %.2e (n=%d).",
				test.statsErrY.mean(), test.statsErrY.std(), test.statsErrY.n() ) );
		System.out.println( String.format( "  - error on Z (~z): %.2e ± %.2e (n=%d).",
				test.statsErrZ.mean(), test.statsErrZ.std(), test.statsErrZ.n() ) );

		test.testRandomAccessBD3rdOrder();
		System.out.println( "\nBackward difference, 3rd order accuracy:" );
		System.out.println( String.format( "  - error on X (~1/x): %.2e ± %.2e (n=%d).",
				test.statsErrX.mean(), test.statsErrX.std(), test.statsErrX.n() ) );
		System.out.println( String.format( "  - error on Y (Cst): %.2e ± %.2e (n=%d).",
				test.statsErrY.mean(), test.statsErrY.std(), test.statsErrY.n() ) );
		System.out.println( String.format( "  - error on Z (~z): %.2e ± %.2e (n=%d).",
				test.statsErrZ.mean(), test.statsErrZ.std(), test.statsErrZ.n() ) );

		test.testRandomAccessBD4thOrder();
		System.out.println( "\nBackward difference, 4th order accuracy:" );
		System.out.println( String.format( "  - error on X (~1/x): %.2e ± %.2e (n=%d).",
				test.statsErrX.mean(), test.statsErrX.std(), test.statsErrX.n() ) );
		System.out.println( String.format( "  - error on Y (Cst): %.2e ± %.2e (n=%d).",
				test.statsErrY.mean(), test.statsErrY.std(), test.statsErrY.n() ) );
		System.out.println( String.format( "  - error on Z (~z): %.2e ± %.2e (n=%d).",
				test.statsErrZ.mean(), test.statsErrZ.std(), test.statsErrZ.n() ) );

		test.testRandomAccessBD5thOrder();
		System.out.println( "\nBackward difference, 5th order accuracy:" );
		System.out.println( String.format( "  - error on X (~1/x): %.2e ± %.2e (n=%d).",
				test.statsErrX.mean(), test.statsErrX.std(), test.statsErrX.n() ) );
		System.out.println( String.format( "  - error on Y (Cst): %.2e ± %.2e (n=%d).",
				test.statsErrY.mean(), test.statsErrY.std(), test.statsErrY.n() ) );
		System.out.println( String.format( "  - error on Z (~z): %.2e ± %.2e (n=%d).",
				test.statsErrZ.mean(), test.statsErrZ.std(), test.statsErrZ.n() ) );

		test.testRandomAccessBD6thOrder();
		System.out.println( "\nBackward difference, 6th order accuracy:" );
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

		// Tests are simple to detect errors, not to benchmark it.
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
	public void testRandomAccessFD1stOrder()
	{
		final RandomAccess< Matrix > ra = GradientRandomAccess.forwardDifference( img, interval, 1 );
		test( ra, 1e-1 );
	}

	@Test
	public void testRandomAccessFD2ndOrder()
	{
		final RandomAccess< Matrix > ra = GradientRandomAccess.forwardDifference( img, interval, 2 );
		test( ra, 1e-2 );
	}

	@Test
	public void testRandomAccessFD3rdOrder()
	{
		final RandomAccess< Matrix > ra = GradientRandomAccess.forwardDifference( img, interval, 3 );
		test( ra, 1e-3 );
	}

	@Test
	public void testRandomAccessFD4thOrder()
	{
		final RandomAccess< Matrix > ra = GradientRandomAccess.forwardDifference( img, interval, 4 );
		test( ra, 1e-4 );
	}

	@Test
	public void testRandomAccessFD5thOrder()
	{
		final RandomAccess< Matrix > ra = GradientRandomAccess.forwardDifference( img, interval, 5 );
		test( ra, 1e-5 );
	}

	@Test
	public void testRandomAccessFD6thOrder()
	{
		final RandomAccess< Matrix > ra = GradientRandomAccess.forwardDifference( img, interval, 6 );
		test( ra, 1e-6 );
	}

	/*
	 * Backward
	 */

	@Test
	public void testRandomAccessBD1stOrder()
	{
		final RandomAccess< Matrix > ra = GradientRandomAccess.backwardDifference( img, interval, 1 );
		test( ra, 1e-1 );
	}

	@Test
	public void testRandomAccessBD2ndOrder()
	{
		final RandomAccess< Matrix > ra = GradientRandomAccess.backwardDifference( img, interval, 2 );
		test( ra, 1e-2 );
	}

	@Test
	public void testRandomAccessBD3rdOrder()
	{
		final RandomAccess< Matrix > ra = GradientRandomAccess.backwardDifference( img, interval, 3 );
		test( ra, 1e-3 );
	}

	@Test
	public void testRandomAccessBD4thOrder()
	{
		final RandomAccess< Matrix > ra = GradientRandomAccess.backwardDifference( img, interval, 4 );
		test( ra, 1e-4 );
	}

	@Test
	public void testRandomAccessBD5thOrder()
	{
		final RandomAccess< Matrix > ra = GradientRandomAccess.backwardDifference( img, interval, 5 );
		test( ra, 1e-5 );
	}

	@Test
	public void testRandomAccessBD6thOrder()
	{
		final RandomAccess< Matrix > ra = GradientRandomAccess.backwardDifference( img, interval, 6 );
		test( ra, 1e-6 );
	}

	/*
	 * Central.
	 */

	@Test
	public void testRandomAccessCD2ndOrder()
	{
		final RandomAccess< Matrix > ra = GradientRandomAccess.centralDifference( img, interval, 2 );
		// Putting a limit on tolerance only makes sense for slowly varying function (from one pixel to the next).
		test( ra, 1e-3 );
	}

	@Test
	public void testRandomAccessCD4thOrder()
	{
		final RandomAccess< Matrix > ra = GradientRandomAccess.centralDifference( img, interval, 4 );
		test( ra, 1e-6 );
	}

	@Test
	public void testRandomAccessCD6thOrder()
	{
		final RandomAccess< Matrix > ra = GradientRandomAccess.centralDifference( img, interval, 6 );
		test( ra, 1e-8 );
	}

	@Test
	public void testRandomAccessCD8thOrder()
	{
		final RandomAccess< Matrix > ra = GradientRandomAccess.centralDifference( img, interval, 8 );
		test( ra, 1e-10 );
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
			final long x = pos[ 0 ];
			final long y = pos[ 1 ];
			final long z = pos[ 2 ];

			ra.setPosition( cursor );
			final Matrix matrix = ra.get();

			ra.localize( posAfter );
			assertArrayEquals( "RandomAccess position changed after access.", pos, posAfter );

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
