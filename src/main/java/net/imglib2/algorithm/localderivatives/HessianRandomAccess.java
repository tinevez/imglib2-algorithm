package net.imglib2.algorithm.localderivatives;

import net.imglib2.Interval;
import net.imglib2.Localizable;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.Sampler;
import net.imglib2.type.numeric.RealType;
import Jama.Matrix;


public class HessianRandomAccess
{


	private HessianRandomAccess()
	{}

	public static final < T extends RealType< T >> RandomAccess< Matrix > centralDifference( final RandomAccessible< T > src, final Interval interval, final int accuracyOrder )
	{
		return new HessianCD2ndOrderRandomAccess< T >( src, interval );
	}


	/*
	 * PRIVATE CLASSES
	 */

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
