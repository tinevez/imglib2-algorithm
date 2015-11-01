package net.imglib2.algorithm.localderivatives;

import net.imglib2.Interval;
import net.imglib2.Localizable;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.Sampler;
import net.imglib2.type.numeric.RealType;
import Jama.Matrix;

/**
 * Static utilities for the Hessian matrix calculation through numerical
 * differentiation.
 * <p>
 * The Hessian matrix is a square matrix containing the second-order partial
 * derivatives of the data sampled by the {@link RandomAccess}. Its elements are
 * such that:
 * <p>
 *
 * <pre>
 * h<sub>i,j</sub> =  &#8706;/&#8706;i &times; &#8706;/&#8706;j
 * </pre>
 * <p>
 * The methods of this class return a {@link RandomAccess} that can be
 * positioned over a {@link RandomAccessible} of type extending {@link RealType}
 * . The {@link RandomAccess#get()} method then returns a
 * <code>n &times; n</code> Jama {@link Matrix}, where <code>n</code> is the
 * dimension of the source.
 *
 * @author Jean-Yves Tinevez - 2015
 */
public class HessianRandomAccess
{
	private HessianRandomAccess()
	{}

	/**
	 * Returns a new {@link RandomAccess}, that can be positioned inside the
	 * specified <code>interval</code>, and can be used to compute the Hessian
	 * matrix by numerical differentiation. The Hessian can be accessed as a
	 * <code>n &times; n</code> Jama {@link Matrix}, where <code>n</code> is the
	 * dimensionality of the source.
	 * <p>
	 * The returned instance is <b>not safe</b> for multithreading: The actual
	 * {@link RandomAccess} returned is used to sample pixel value when
	 * calculating the gradient. The Jama matrix returned by the
	 * {@link RandomAccess#get()} method is the same at each call.
	 * <p>
	 * The gradient calculation involves numerical differentiation by central
	 * differences. Supported accuracy orders are 2 and 4. The actual interval
	 * effectively sampled by the gradient calculation is expanded by the
	 * central differences differentiation. Borders of 1 and 2 pixels are added
	 * for accuracy orders 2 and 4 respectively on all sides of all dimensions.
	 * It is the responsibility of the caller to ensure that the source is
	 * accessible over this expanded interval.
	 *
	 * @param src
	 *            the source.
	 * @param interval
	 *            the interval on which the returned {@link RandomAccess} will
	 *            be positioned.
	 * @param accuracyOrder
	 *            the desired accuracy orders. Orders 2 and 4 are supported.
	 * @return a new {@link RandomAccess}.
	 * @throws IllegalArgumentException
	 *             if the specified accuracy order is not 2 or 4.
	 */
	public static final < T extends RealType< T >> RandomAccess< Matrix > centralDifference( final RandomAccessible< T > src, final Interval interval, final int accuracyOrder )
	{
		switch ( accuracyOrder )
		{
		case 2:
			return new HessianCD2ndOrderRandomAccess< T >( src, interval );
		case 4:
			return new HessianCD4thOrderRandomAccess< T >( src, interval );
		default:
			throw new IllegalArgumentException( "Accuracy order can only be 2 or 4. Was: " + accuracyOrder + "." );
		}
	}

	/**
	 * Returns a new {@link RandomAccess}, that can be positioned inside the
	 * specified <code>interval</code>, and can be used to compute the Hessian
	 * matrix by numerical differentiation. The Hessian can be accessed as a
	 * <code>n &times; n</code> Jama {@link Matrix}, where <code>n</code> is the
	 * dimensionality of the source.
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

	/*
	 * PRIVATE CLASSES
	 */

	private static final class HessianCD4thOrderRandomAccess< T extends RealType< T >> extends AbstractHessianRandomAccess< T >
	{

		private HessianCD4thOrderRandomAccess( final RandomAccessible< T > src, final Interval interval )
		{
			super( src, interval );
		}

		@Override
		public Matrix get()
		{
			final double a2 = randomAccess.get().getRealDouble();
			for ( int d = 0; d < n; ++d )
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

				matrix.set( d, d, -1. / 12. * ( a4 + a0 ) + 4. / 3. * ( a3 + a1 ) - 5. / 2. * a2 );

				for ( int e = d + 1; e < n; ++e )
				{
					// We start from center point.
					randomAccess.move( -2, d );
					randomAccess.move( -2, e );
					final double a0b0 = randomAccess.get().getRealDouble();
					randomAccess.fwd( d );
					final double a1b0 = randomAccess.get().getRealDouble();
					randomAccess.move( 2, d );
					final double a3b0 = randomAccess.get().getRealDouble();
					randomAccess.fwd( d );
					final double a4b0 = randomAccess.get().getRealDouble();

					randomAccess.fwd( e );
					final double a4b1 = randomAccess.get().getRealDouble();
					randomAccess.bck( d );
					final double a3b1 = randomAccess.get().getRealDouble();
					randomAccess.move( -2, d );
					final double a1b1 = randomAccess.get().getRealDouble();
					randomAccess.bck( d );
					final double a0b1 = randomAccess.get().getRealDouble();

					randomAccess.move( 2, e );
					final double a0b3 = randomAccess.get().getRealDouble();
					randomAccess.fwd( d );
					final double a1b3 = randomAccess.get().getRealDouble();
					randomAccess.move( 2, d );
					final double a3b3 = randomAccess.get().getRealDouble();
					randomAccess.fwd( d );
					final double a4b3 = randomAccess.get().getRealDouble();

					randomAccess.fwd( e );
					final double a4b4 = randomAccess.get().getRealDouble();
					randomAccess.bck( d );
					final double a3b4 = randomAccess.get().getRealDouble();
					randomAccess.move( -2, d );
					final double a1b4 = randomAccess.get().getRealDouble();
					randomAccess.bck( d );
					final double a0b4 = randomAccess.get().getRealDouble();

					final double v =
							1. / 12. * ( -1. / 12. * ( a0b4 - a0b0 ) + 2. / 3. * ( a0b3 - a0b1 ) )
									- 1. / 12. * ( -1. / 12. * ( a4b4 - a4b0 ) + 2. / 3. * ( a4b3 - a4b1 ) )
									+ 2. / 3. * ( -1. / 12. * ( a3b4 - a3b0 ) + 2. / 3. * ( a3b3 - a3b1 ) )
									- 2. / 3. * ( -1. / 12. * ( a1b4 - a1b0 ) + 2. / 3. * ( a1b3 - a1b1 ) );
					matrix.set( d, e, v );
					matrix.set( e, d, v );

					// Move back to the center.
					randomAccess.move( 2, d );
					randomAccess.move( -2, e );
				}
			}

			return matrix;
		}

		@Override
		public RandomAccess< Matrix > copyRandomAccess()
		{
			return new HessianCD4thOrderRandomAccess< T >( src, interval );
		}

		@Override
		public Sampler< Matrix > copy()
		{
			return copyRandomAccess();
		}
	}

	private static final class HessianCD2ndOrderRandomAccess< T extends RealType< T >> extends AbstractHessianRandomAccess< T >
	{

		private HessianCD2ndOrderRandomAccess( final RandomAccessible< T > src, final Interval interval )
		{
			super( src, interval );
		}

		@Override
		public Matrix get()
		{
			final double a1 = randomAccess.get().getRealDouble();
			for ( int d = 0; d < n; ++d )
			{
				randomAccess.bck( d );
				final double a0 = randomAccess.get().getRealDouble();
				randomAccess.move( 2, d );
				final double a2 = randomAccess.get().getRealDouble();
				matrix.set( d, d, a2 - 2 * a1 + a0 );
				randomAccess.bck( d );

				for ( int e = d + 1; e < n; ++e )
				{
					// We start from center point.
					randomAccess.fwd( d );
					randomAccess.fwd( e );
					final double a2b2 = randomAccess.get().getRealDouble();
					randomAccess.move( -2, d );
					final double a0b2 = randomAccess.get().getRealDouble();
					randomAccess.move( -2, e );
					final double a0b0 = randomAccess.get().getRealDouble();
					randomAccess.move( 2, d );
					final double a2b0 = randomAccess.get().getRealDouble();
					// back to the original position
					randomAccess.bck( d );
					randomAccess.fwd( e );
					final double v = ( a2b2 - a0b2 - a2b0 + a0b0 ) * 0.25;
					matrix.set( d, e, v );
					matrix.set( e, d, v );
				}
			}

			return matrix;
		}

		@Override
		public RandomAccess< Matrix > copyRandomAccess()
		{
			return new HessianCD2ndOrderRandomAccess< T >( src, interval );
		}

		@Override
		public Sampler< Matrix > copy()
		{
			return copyRandomAccess();
		}

	}

	private static abstract class AbstractHessianRandomAccess< T extends RealType< T >> implements RandomAccess< Matrix >
	{
		protected final RandomAccess< T > randomAccess;

		protected final Matrix matrix;

		protected final RandomAccessible< T > src;

		protected final Interval interval;

		protected final int n;

		protected AbstractHessianRandomAccess( final RandomAccessible< T > src, final Interval interval )
		{
			this.src = src;
			this.interval = interval;
			this.n = src.numDimensions();
			randomAccess = src.randomAccess( interval );
			matrix = new Matrix( n, n );
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
