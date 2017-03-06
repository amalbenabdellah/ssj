package umontreal.ssj.stat.density;

import umontreal.ssj.stat.KernelDensityEstimator1d;

/**
 * Density estimator based on Gaussian kernel, as proposed by Botev.
 * 
 * @author Lecuyer
 *
 */
public class DEGaussianKDEBotev implements DensityEstimator {
	double a,b,bandwidthKDE;
	int m;
	KernelDensityEstimator1d kde;

	DEGaussianKDEBotev (double a, double b) {
		this.a=a;
		this.b=b;
		
		
	}
    DEGaussianKDEBotev (double a, double b,int m, double bandwidthKDE ) {
    	this.a=a;
		this.b=b;
		this.m=m;
		this.bandwidthKDE=bandwidthKDE;
	}

	
	public void setRange(double a, double b) {
		// TODO Auto-generated method stub
		this.a=a;
		this.b=b;
	}

	
	public void constructDensity(double[] x) {
		// TODO Auto-generated method stub
		 kde =new KernelDensityEstimator1d ();
		 kde.kde(x,m,a,b,bandwidthKDE);
	}

	
	public double evalDensity(double x) {
		// TODO Auto-generated method stub
		double[] densi=kde.getDensity();
		double h=(b-a)/m;
		return densi[(int) ((x -a) /h)];
	}

	
	public void evalDensity(double[] x, double[] f) {
		// TODO Auto-generated method stub
		for(int i=0; i < x.length; i++)
			f[i]=evalDensity(x[i]);
	}
		
}
