package net.imglib2.algorithm.localderivatives;

class MeanStd
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