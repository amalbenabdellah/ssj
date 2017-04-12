package umontreal.ssj.stat.density;

import umontreal.ssj.probdist.NormalDist;
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

	public DEGaussianKDEBotev (double a, double b) {
		
		this.a=a;
		this.b=b;
		
		
	}
	
    public DEGaussianKDEBotev (double a, double b, int m, double bandwidthKDE ) {
    	
    	this.a = a;
		this.b = b;
		this.m = m;
		this.bandwidthKDE = bandwidthKDE;
		kde = null;
		
	}
    
    public void setRange(double a, double b) {
		
		this.a = a;
		this.b = b;
	}

	
	public double getA() {
		return a;
	}
	public void setA(double a) {
		this.a = a;
	}
	public double getB() {
		return b;
	}
	public void setB(double b) {
		this.b = b;
	}
	public double getBandwidthKDE() {
		return bandwidthKDE;
	}
	public void setBandwidthKDE(double bandwidthKDE) {
		this.bandwidthKDE = bandwidthKDE;
	}
	public int getM() {
		return m;
	}
	public void setM(int m) {
		this.m = m;
	}
	public KernelDensityEstimator1d getKde() {
		return kde;
	}
	public void setKde(KernelDensityEstimator1d kde) {
		this.kde = kde;
	}
	

	
	// Estimate the density
	
	public void constructDensity(double[] x) {
		
		 kde = new KernelDensityEstimator1d ();
		 kde.kde(x,m,a,b,bandwidthKDE);
		
	}

	public double evalDensity(double x, int n) {
		return 0.0;
	}
	// Evaluate  the density  at a point
	
	public double evalDensity(double x, double[] data) {
		
		double sum=0.0;
		for(int i=0; i < data.length; i++){
		sum += NormalDist.density(0, 1, (data[i]-x)/bandwidthKDE);
		}
		return sum/(bandwidthKDE*data.length);
		
	/*	double[] densi = kde.getDensity();
		double h = (b-a)/m;
		return densi[(int) ((x -a) /h)];*/
		
	}

	
	
	
    /*public void evalDensity(double[] x, double[] f) {
		
		for(int i=0; i < x.length; i++)
			f[i] = evalDensity(x[i],x.length);
	}*/
	
	// Evaluate the density at eval points
	
	public void evalDensity(double[] x, double[] f,double[] data) {
		
		for(int i=0; i < x.length; i++)
			f[i] = evalDensity(x[i], data);
	}
		
}
