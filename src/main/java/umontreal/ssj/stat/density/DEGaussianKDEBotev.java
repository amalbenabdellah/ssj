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
	double a,b,h;
	int m;
	KernelDensityEstimator1d kde;

	public DEGaussianKDEBotev (double a, double b) {
		
		this.a=a;
		this.b=b;
		
		
	}
	
    public DEGaussianKDEBotev (double a, double b, int m, double h ) {
    	
    	this.a = a;
		this.b = b;
		this.m = m;
		this.h = h;
		kde = null;
		
	}
    
    public DEGaussianKDEBotev (double a, double b, int m ) {
    	
    	this.a = a;
		this.b = b;
		this.m = m;
		kde = null;
		
	}
    
    public void setRange(double a, double b) {
		
		this.a = a;
		this.b = b;
	}

	
    public void seth(double h) {
		this.h = h;
	}
	public double getA() {
		return a;
	}
	public double getB() {
		return b;
	}
	
	
	public int getM() {
		return m;
	}
	public double geth(){
		return h;
	}
	
	
	

	
	/**	
     * Constructs a density estimator from 
     * the data points in vector x.
     * @param x  data points.
     */
	
	public void constructDensity(double[] x) {
		
		 kde = new KernelDensityEstimator1d ();
		 
		 //For fixed bandwidth
		// kde.kde(x,m,a,b,bandwidthKDE);
		//For variable bandwidth
		 kde.kde(x,m,a,b);
	}

	/**
	 * Returns the value of the density evaluated at x.
	 */
	public double evalDensity(double x, int n) {
		double[] densi = kde.getDensity();
		double h = (b-a)/m;
		return densi[(int) ((x -a) /h)];
	}
	/**
	 * Returns the value of the density evaluated at x.
	 */
	public double evalDensity(double x, double[] data) {
		
		
		double[] densi = kde.getDensity();
		double h = (b-a)/m;
		return densi[(int) ((x -a) /h)];
		
	}

	
	/**	
    * Returns in array density the value of the density at the evaluation points in x.
    * These two arrays must have the same size. 
    * @param x  evaluation points
    * @param density  values of the density
    */
	public void evalDensity(double[] x, double[] f,double[] data) {
		
		for(int i=0; i < x.length; i++)
			f[i] = evalDensity(x[i], data);
	}
	
	public String toString(){
		return "Kernel Botev _"+ geth();
	}
		
}
