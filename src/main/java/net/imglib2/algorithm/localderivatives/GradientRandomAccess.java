package net.imglib2.algorithm.localderivatives;

import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.Localizable;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.type.numeric.RealType;
import net.imglib2.util.Intervals;
import Jama.Matrix;

/**
 * Static utilities for gradient vector calculation through numerical
 * differentiation.
 * <p>
 * The methods of this class return a {@link RandomAccess} that can be
 * positioned over a {@link RandomAccessible} of type extending {@link RealType}
 * . The {@link RandomAccess#get()} method then returns a
 * <code>n &times; 1</code>Jama {@link Matrix}, where <code>n</code> is the
 * dimension of the source.
 *
 * @author Jean-Yves Tinevez - 2015
 *
 */
public class GradientRandomAccess
{
	private GradientRandomAccess()
	{}

	/**
	 * Returns a new {@link RandomAccess}, that can be positioned inside the
	 * specified <code>interval</code>, and can be used to compute the local
	 * gradient by numerical differentiation.
	 * <p>
	 * The returned instance is <b>not safe</b> for multithreading: The actual
	 * {@link RandomAccess} returned is used to sample pixel value when
	 * calculating the gradient. The Jama matrix returned by the
	 * {@link RandomAccess#get()} method is the same at each call.
	 * <p>
	 * The gradient calculation involves numerical differentiation by central
	 * differences with an accuracy order of 2. The actual interval effectively
	 * sampled by the gradient calculation is expanded by the central
	 * differences differentiation. Borders of 1 pixel are added on all sides of
	 * all dimensions. It is the responsibility of the caller to ensure that the
	 * source is accessible over this expanded interval.
	 *
	 * @param src
	 *            the source.
	 * @param interval
	 *            the interval on which the returned {@link RandomAccess} will
	 *            be positioned.
	 * @return a new {@link RandomAccess}.
	 */
	public static final < T extends RealType< T >> RandomAccess< Matrix > centralDifference( final RandomAccessible< T > src, final Interval interval )
	{
		return centralDifference( src, interval, 2 );
	}

	/**
	 * Returns a new {@link RandomAccess}, that can be positioned inside the
	 * specified <code>interval</code>, and can be used to compute the local
	 * gradient by numerical differentiation.
	 * <p>
	 * The returned instance is <b>not safe</b> for multithreading: The actual
	 * {@link RandomAccess} returned is used to sample pixel value when
	 * calculating the gradient. The Jama matrix returned by the
	 * {@link RandomAccess#get()} method is the same at each call.
	 * <p>
	 * The gradient calculation involves numerical differentiation by central
	 * differences. Supported accuracy orders are 2, 4,6 and 8. The actual
	 * interval effectively sampled by the gradient calculation is expanded by
	 * the central differences differentiation. Borders of 1, 2, 3 and 4 pixels
	 * are added for accuracy orders 2, 4, 6 and 8 respectively on all sides of
	 * all dimensions. It is the responsibility of the caller to ensure that the
	 * source is accessible over this expanded interval.
	 *
	 * @param src
	 *            the source.
	 * @param interval
	 *            the interval on which the returned {@link RandomAccess} will
	 *            be positioned.
	 * @param accuracyOrder
	 *            the desired accuracy orders. Orders 2, 4, 6 and 8 are
	 *            supported.
	 * @return a new {@link RandomAccess}.
	 * @throws IllegalArgumentException
	 *             if the specified accuracy order is not 2, 4, 6 or 8.
	 */
	public static final < T extends RealType< T >> RandomAccess< Matrix > centralDifference( final RandomAccessible< T > src, final Interval interval, final int accuracyOrder )
	{
		final FinalInterval expandedInterval = Intervals.expand( interval, accuracyOrder / 2 );
		switch ( accuracyOrder )
		{
		case 2:
			return new GradientCD2ndOrderRandomAccess< T >( src, expandedInterval );
		case 4:
			return new GradientCD4thOrderRandomAccess< T >( src, expandedInterval );
		case 6:
			return new GradientCD6thOrderRandomAccess< T >( src, expandedInterval );
		case 8:
			return new GradientCD8thOrderRandomAccess< T >( src, expandedInterval );
		default:
			throw new IllegalArgumentException( "Accuracy order can only be 2, 4, 6 or 8. Was: " + accuracyOrder + "." );
		}
	}

	public static final < T extends RealType< T >> RandomAccess< Matrix > forwardDifference( final RandomAccessible< T > src, final Interval interval )
	{
		return forwardDifference( src, interval, 2 );
	}

	public static final < T extends RealType< T >> RandomAccess< Matrix > forwardDifference( final RandomAccessible< T > src, final Interval interval, final int accuracyOrder )
	{
		final int n = src.numDimensions();
		final long[] minmax = new long[ 2 * n ];
		for ( int d = 0; d < n; d++ )
		{
			minmax[ d ] = interval.min( d );
			minmax[ n + d ] = interval.max( d ) + accuracyOrder;
		}
		final FinalInterval expandedInterval = Intervals.createMinMax( minmax );
		switch ( accuracyOrder )
		{
		case 1:
			return new GradientFD1stOrderRandomAccess< T >( src, expandedInterval );
		case 2:
			return new GradientFD2nOrderRandomAccess< T >( src, expandedInterval );
		case 3:
			return new GradientFD3rdOrderRandomAccess< T >( src, expandedInterval );
		case 4:
			return new GradientFD4thOrderRandomAccess< T >( src, expandedInterval );
		case 5:
			return new GradientFD5thOrderRandomAccess< T >( src, expandedInterval );
		case 6:
			return new GradientFD6thOrderRandomAccess< T >( src, expandedInterval );
		default:
			throw new IllegalArgumentException( "Accuracy order can only be 1, 2, 3, 4, 5 or 6. Was: " + accuracyOrder + "." );
		}
	}

	private static class GradientFD6thOrderRandomAccess< T extends RealType< T >> extends AbstractGradientRandomAccess< T >
	{

		private GradientFD6thOrderRandomAccess( final RandomAccessible< T > src, final Interval interval )
		{
			super( src, interval );
		}

		@Override
		public Matrix get()
		{
			for ( int d = 0; d < numDimensions(); d++ )
			{
				randomAccess.move( 6, d );
				final double a6 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				final double a5 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				final double a4 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				final double a3 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				final double a2 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				final double a1 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				final double a0 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				matrix.set( d, 0, -1. / 6. * a6 + 6. / 5. * a5 - 15. / 4. * a4 + 20. / 3. * a3 - 15. / 2. * a2 + 6. * a1 - 49. / 20. * a0 );
			}
			return matrix;
		}

		@Override
		public RandomAccess< Matrix > copyRandomAccess()
		{
			return new GradientFD6thOrderRandomAccess< T >( src, interval );
		}

		@Override
		public RandomAccess< Matrix > copy()
		{
			return copyRandomAccess();
		}
	}

	private static class GradientFD5thOrderRandomAccess< T extends RealType< T >> extends AbstractGradientRandomAccess< T >
	{

		private GradientFD5thOrderRandomAccess( final RandomAccessible< T > src, final Interval interval )
		{
			super( src, interval );
		}

		@Override
		public Matrix get()
		{
			for ( int d = 0; d < numDimensions(); d++ )
			{
				randomAccess.move( 5, d );
				final double a5 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				final double a4 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				final double a3 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				final double a2 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				final double a1 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				final double a0 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				matrix.set( d, 0, 1. / 5. * a5 - 5. / 4. * a4 + 10. / 3. * a3 - 5. * a2 + 5. * a1 - 137. / 60. * a0 );
			}
			return matrix;
		}

		@Override
		public RandomAccess< Matrix > copyRandomAccess()
		{
			return new GradientFD5thOrderRandomAccess< T >( src, interval );
		}

		@Override
		public RandomAccess< Matrix > copy()
		{
			return copyRandomAccess();
		}
	}

	private static class GradientFD4thOrderRandomAccess< T extends RealType< T >> extends AbstractGradientRandomAccess< T >
	{

		private GradientFD4thOrderRandomAccess( final RandomAccessible< T > src, final Interval interval )
		{
			super( src, interval );
		}

		@Override
		public Matrix get()
		{
			for ( int d = 0; d < numDimensions(); d++ )
			{
				randomAccess.move( 4, d );
				final double a4 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				final double a3 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				final double a2 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				final double a1 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				final double a0 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				matrix.set( d, 0, -1. / 4. * a4 + 4. / 3. * a3 - 3. * a2 + 4. * a1 - 25. / 12. * a0 );
			}
			return matrix;
		}

		@Override
		public RandomAccess< Matrix > copyRandomAccess()
		{
			return new GradientFD4thOrderRandomAccess< T >( src, interval );
		}

		@Override
		public RandomAccess< Matrix > copy()
		{
			return copyRandomAccess();
		}
	}

	private static class GradientFD3rdOrderRandomAccess< T extends RealType< T >> extends AbstractGradientRandomAccess< T >
	{

		private GradientFD3rdOrderRandomAccess( final RandomAccessible< T > src, final Interval interval )
		{
			super( src, interval );
		}

		@Override
		public Matrix get()
		{
			for ( int d = 0; d < numDimensions(); d++ )
			{
				randomAccess.move( 3, d );
				final double a3 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				final double a2 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				final double a1 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				final double a0 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				matrix.set( d, 0, 1. / 3. * a3 - 3. / 2. * a2 + 3. * a1 - 11. / 6. * a0 );
			}
			return matrix;
		}

		@Override
		public RandomAccess< Matrix > copyRandomAccess()
		{
			return new GradientFD3rdOrderRandomAccess< T >( src, interval );
		}

		@Override
		public RandomAccess< Matrix > copy()
		{
			return copyRandomAccess();
		}
	}

	private static class GradientFD2nOrderRandomAccess< T extends RealType< T >> extends AbstractGradientRandomAccess< T >
	{

		private GradientFD2nOrderRandomAccess( final RandomAccessible< T > src, final Interval interval )
		{
			super( src, interval );
		}

		@Override
		public Matrix get()
		{
			for ( int d = 0; d < numDimensions(); d++ )
			{
				randomAccess.move( 2, d );
				final double a2 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				final double a1 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				final double a0 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				matrix.set( d, 0, -1. / 2. * a2 + 2. * a1 - 3. / 2. * a0 );
			}
			return matrix;
		}

		@Override
		public RandomAccess< Matrix > copyRandomAccess()
		{
			return new GradientFD2nOrderRandomAccess< T >( src, interval );
		}

		@Override
		public RandomAccess< Matrix > copy()
		{
			return copyRandomAccess();
		}
	}

	private static class GradientFD1stOrderRandomAccess< T extends RealType< T >> extends AbstractGradientRandomAccess< T >
	{

		private GradientFD1stOrderRandomAccess( final RandomAccessible< T > src, final Interval interval )
		{
			super( src, interval );
		}

		@Override
		public Matrix get()
		{
			for ( int d = 0; d < numDimensions(); d++ )
			{
				final double a0 = randomAccess.get().getRealDouble();
				randomAccess.fwd( d );
				final double a1 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				matrix.set( d, 0, a1 - a0 );
			}
			return matrix;
		}

		@Override
		public RandomAccess< Matrix > copyRandomAccess()
		{
			return new GradientFD1stOrderRandomAccess< T >( src, interval );
		}

		@Override
		public RandomAccess< Matrix > copy()
		{
			return copyRandomAccess();
		}
	}

	private static class GradientCD2ndOrderRandomAccess< T extends RealType< T >> extends AbstractGradientRandomAccess< T >
	{

		private GradientCD2ndOrderRandomAccess( final RandomAccessible< T > src, final Interval interval )
		{
			super( src, interval );
		}

		@Override
		public Matrix get()
		{
			for ( int d = 0; d < numDimensions(); d++ )
			{
				randomAccess.bck( d );
				final double a0 = randomAccess.get().getRealDouble();
				randomAccess.move( 2, d );
				final double a2 = randomAccess.get().getRealDouble();
				randomAccess.bck( d );
				matrix.set( d, 0, 1. / 2. * ( a2 - a0 ) );
			}
			return matrix;
		}

		@Override
		public RandomAccess< Matrix > copyRandomAccess()
		{
			return new GradientCD2ndOrderRandomAccess< T >( src, interval );
		}

		@Override
		public RandomAccess< Matrix > copy()
		{
			return copyRandomAccess();
		}
	}

	private static class GradientCD4thOrderRandomAccess< T extends RealType< T >> extends AbstractGradientRandomAccess< T >
	{

		private GradientCD4thOrderRandomAccess( final RandomAccessible< T > src, final Interval interval )
		{
			super( src, interval );
		}

		@Override
		public Matrix get()
		{
			for ( int d = 0; d < numDimensions(); d++ )
			{
				randomAccess.move( -2, d );
				final double a0 = randomAccess.get().getRealDouble();
				randomAccess.fwd( d );
				final double a1 = randomAccess.get().getRealDouble();
				randomAccess.move( 2, d );
				final double a3 = randomAccess.get().getRealDouble();
				randomAccess.fwd( d );
				final double a4 = randomAccess.get().getRealDouble();
				randomAccess.move( -2, d );
				matrix.set( d, 0, -1. / 12. * ( a4 - a0 ) + 2. / 3. * ( a3 - a1 ) );
			}
			return matrix;
		}

		@Override
		public RandomAccess< Matrix > copyRandomAccess()
		{
			return new GradientCD4thOrderRandomAccess< T >( src, interval );
		}

		@Override
		public RandomAccess< Matrix > copy()
		{
			return copyRandomAccess();
		}
	}

	private static class GradientCD6thOrderRandomAccess< T extends RealType< T >> extends AbstractGradientRandomAccess< T >
	{

		private GradientCD6thOrderRandomAccess( final RandomAccessible< T > src, final Interval interval )
		{
			super( src, interval );
		}

		@Override
		public Matrix get()
		{
			for ( int d = 0; d < numDimensions(); d++ )
			{
				randomAccess.move( -3, d );
				final double a0 = randomAccess.get().getRealDouble();
				randomAccess.fwd( d );
				final double a1 = randomAccess.get().getRealDouble();
				randomAccess.fwd( d );
				final double a2 = randomAccess.get().getRealDouble();
				randomAccess.move( 2, d );
				final double a4 = randomAccess.get().getRealDouble();
				randomAccess.fwd( d );
				final double a5 = randomAccess.get().getRealDouble();
				randomAccess.fwd( d );
				final double a6 = randomAccess.get().getRealDouble();
				randomAccess.move( -3, d );
				matrix.set( d, 0, 1. / 60. * ( a6 - a0 ) - 3. / 20. * ( a5 - a1 ) + 3. / 4. * ( a4 - a2 ) );
			}
			return matrix;
		}

		@Override
		public RandomAccess< Matrix > copyRandomAccess()
		{
			return new GradientCD6thOrderRandomAccess< T >( src, interval );
		}

		@Override
		public RandomAccess< Matrix > copy()
		{
			return copyRandomAccess();
		}
	}

	private static class GradientCD8thOrderRandomAccess< T extends RealType< T >> extends AbstractGradientRandomAccess< T >
	{

		private GradientCD8thOrderRandomAccess( final RandomAccessible< T > src, final Interval interval )
		{
			super( src, interval );
		}

		@Override
		public Matrix get()
		{
			for ( int d = 0; d < numDimensions(); d++ )
			{
				randomAccess.move( -4, d );
				final double a0 = randomAccess.get().getRealDouble();
				randomAccess.fwd( d );
				final double a1 = randomAccess.get().getRealDouble();
				randomAccess.fwd( d );
				final double a2 = randomAccess.get().getRealDouble();
				randomAccess.fwd( d );
				final double a3 = randomAccess.get().getRealDouble();
				randomAccess.move( 2, d );
				final double a5 = randomAccess.get().getRealDouble();
				randomAccess.fwd( d );
				final double a6 = randomAccess.get().getRealDouble();
				randomAccess.fwd( d );
				final double a7 = randomAccess.get().getRealDouble();
				randomAccess.fwd( d );
				final double a8 = randomAccess.get().getRealDouble();
				randomAccess.move( -4, d );
				matrix.set( d, 0, -1. / 280. * ( a8 - a0 ) + 4. / 105. * ( a7 - a1 ) - 1. / 5. * ( a6 - a2 ) + 4. / 5. * ( a5 - a3 ) );
			}
			return matrix;
		}

		@Override
		public RandomAccess< Matrix > copyRandomAccess()
		{
			return new GradientCD8thOrderRandomAccess< T >( src, interval );
		}

		@Override
		public RandomAccess< Matrix > copy()
		{
			return copyRandomAccess();
		}
	}

	private static abstract class AbstractGradientRandomAccess< T extends RealType< T >> implements RandomAccess< Matrix >
	{
		protected final RandomAccess< T > randomAccess;

		protected final Matrix matrix;

		protected final RandomAccessible< T > src;

		protected final Interval interval;

		protected AbstractGradientRandomAccess( final RandomAccessible< T > src, final Interval interval )
		{
			this.src = src;
			this.interval = interval;
			randomAccess = src.randomAccess( interval );
			matrix = new Matrix( src.numDimensions(), 1 );
		}

		@Override
		public void localize( final int[] position )
		{
			randomAccess.localize( position );
		}

		@Override
		public void localize( final long[] position )
		{
			randomAccess.localize( position );
		}

		@Override
		public void localize( final float[] position )
		{
			randomAccess.localize( position );
		}

		@Override
		public void localize( final double[] position )
		{
			randomAccess.localize( position );
		}

		@Override
		public int getIntPosition( final int d )
		{
			return randomAccess.getIntPosition( d );
		}

		@Override
		public long getLongPosition( final int d )
		{
			return randomAccess.getLongPosition( d );
		}

		@Override
		public float getFloatPosition( final int d )
		{
			return randomAccess.getFloatPosition( d );
		}

		@Override
		public double getDoublePosition( final int d )
		{
			return randomAccess.getDoublePosition( d );
		}

		@Override
		public int numDimensions()
		{
			return randomAccess.numDimensions();
		}

		@Override
		public void fwd( final int d )
		{
			randomAccess.fwd( d );
		}

		@Override
		public void bck( final int d )
		{
			randomAccess.bck( d );
		}

		@Override
		public void move( final int distance, final int d )
		{
			randomAccess.move( distance, d );
		}

		@Override
		public void move( final long distance, final int d )
		{
			randomAccess.move( distance, d );
		}

		@Override
		public void move( final Localizable localizable )
		{
			randomAccess.move( localizable );
		}

		@Override
		public void move( final int[] distance )
		{
			randomAccess.move( distance );
		}

		@Override
		public void move( final long[] distance )
		{
			randomAccess.move( distance );
		}

		@Override
		public void setPosition( final Localizable localizable )
		{
			randomAccess.setPosition( localizable );
		}

		@Override
		public void setPosition( final int[] position )
		{
			randomAccess.setPosition( position );
		}

		@Override
		public void setPosition( final long[] position )
		{
			randomAccess.setPosition( position );
		}

		@Override
		public void setPosition( final int position, final int d )
		{
			randomAccess.setPosition( position, d );
		}

		@Override
		public void setPosition( final long position, final int d )
		{
			randomAccess.setPosition( position, d );
		}
	}
}
