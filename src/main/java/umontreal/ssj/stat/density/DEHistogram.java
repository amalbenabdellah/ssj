package umontreal.ssj.stat.density;
import umontreal.ssj.stat.ScaledHistogram;
import umontreal.ssj.stat.TallyHistogram;

/**
* This class provides methods to compute a histogram density estimator from a
* set of @f$n@f$ individual observations @f$x_0, …, x_{n-1}@f$, and returns
* its value at a  set of selected points. 
*/


public class DEHistogram implements DensityEstimator {
	double a,b, h;
	int m;
	ScaledHistogram histDensity;

	
	
	
	public DEHistogram (double a, double b) {
		this.a = a;
		this.b = b;

		histDensity = null;
		
	}
    public DEHistogram (double a, double b,int  m) {
    	this.a = a;
		this.b = b;
		this.m = m;
		this.h= b-a/m;
		histDensity = null;
	}

	
	public void setRange(double a, double b) {
		
		this.a = a;
		this.b = b;
	}
	
	
	
	
	
	
	public void seth(double h) {
		this.m = (int) ((b-a)/h);
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
	
	public double geth() {
		return h;
	}
	
	
	

	 /**	
     * Constructs a density estimator from 
     * the data points in vector x.
     * @param x  data points.
     */
	
	public void constructDensity(double[] x) {
		
		TallyHistogram hist = new TallyHistogram(a,b,m);
		hist.fillFromArray(x);
		histDensity = new ScaledHistogram(hist,1.0);
		
	}

	/**
	 * Returns the value of the density evaluated at x.
	 */
	
	public double evalDensity(double x, int n ) {
		double h=(b-a)/m;
		
		return histDensity.getHeight(x)/(n*h);
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
		f[i] = evalDensity(x[i],data.length);
	}
	
	
	public String toString(){
		return "Histogram with Numbins"+ getM();
	}
	
	
}
