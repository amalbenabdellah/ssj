package umontreal.ssj.stat.density;

import umontreal.ssj.stat.ScaledHistogram;
import umontreal.ssj.stat.TallyHistogram;

/**
* This class provides methods to compute an averaged shifted histogram with linear weight density estimator from a
* set of @f$n@f$ individual observations @f$x_0, …, x_{n-1}@f$, and returns
* its value at a  set of selected points. 
*/
public class  DEAveragedShiftedHistogram  implements DensityEstimator {
	double a,b;
	int m, r;
	ScaledHistogram histDensity;
	ScaledHistogram histAsh;

	
	
	public DEAveragedShiftedHistogram (double a, double b) {
		this.a = a;
		this.b = b;
		histDensity = null;
		histAsh = null;
		
		
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
	public void seth(double h) {
		this.m = (int) ((b-a)/h);
	}
	public double geth() {
		return (b-a)/m;
	}
	
	public DEAveragedShiftedHistogram (double a, double b, int m, int r) {
	this.a = a;
	this.b = b;
	this.m = m;
	this.r = r;
	histDensity = null;
	histAsh = null;

	}
	public void setRange(double a, double b) {
		
		this.a = a;
		this.b = b;
		histDensity = null;
		histAsh = null;
	}

	
	 /**	
     * Constructs a density estimator from 
     * the data points in vector x.
     * @param x  data points.
     */
	public void constructDensity(double[] x) {
		
		TallyHistogram hist = new TallyHistogram(a,b,m);
		hist.fillFromArray(x);
		
		histDensity=new ScaledHistogram(hist,1.0);
	    histAsh=histDensity.averageShiftedHistogramTrunc(r);
		histAsh.rescale(1.0);
	}

	/**
	 * Returns the value of the density evaluated at x.
	 */
	
	public double evalDensity(double x, int n) {
		double h = (b - a)/ r*m;	
		double evaldens = histAsh.getHeights()[(int) ((x - a)/h)]/(n*h*r*r);
		return evaldens;
	}
	
 
	/**	
     * Returns in array density the value of the density at the evaluation points in x.
     * These two arrays must have the same size. 
     * @param x  evaluation points
     * @param density  values of the density
     * @param data array of observations
     */
	
    public void evalDensity(double[] x, double[] f, double[] data) {
    	
    	for(int i=0; i < x.length; i++)
			f[i] = evalDensity(x[i], data.length);
	}
    
    public String toString(){
		return "ASH linear_"+ getM();
	}
	
}
