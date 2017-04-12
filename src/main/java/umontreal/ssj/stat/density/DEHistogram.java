package umontreal.ssj.stat.density;
import umontreal.ssj.stat.ScaledHistogram;
import umontreal.ssj.stat.TallyHistogram;

/**
 * 
 * @author Lecuyer
 *
 * A histogram used a a density estimator.
 * 
 */
public class DEHistogram implements DensityEstimator {
	double a,b;
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
		histDensity = null;
	}

	
	public void setRange(double a, double b) {
		
		this.a = a;
		this.b = b;
		histDensity = null;
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
	public int getM() {
		return m;
	}
	public void setM(int m) {
		this.m = m;
	}
	public ScaledHistogram getHistDensity() {
		return histDensity;
	}
	public void setHistDensity(ScaledHistogram histDensity) {
		this.histDensity = histDensity;
	}

	// Estimate the density
	
	public void constructDensity(double[] x) {
		
		TallyHistogram hist = new TallyHistogram(a,b,m);
		hist.fillFromArray(x);
		histDensity = new ScaledHistogram(hist,1.0);
		
	}

	// Evaluate  the density  at a point
	
	public double evalDensity(double x, int n ) {
		double h=(b-a)/m;
		
		return histDensity.getHeight(x)/(n*h);
	}

	// Evaluate the density at eval points
	public void evalDensity(double[] x, double[] f, double[] data) {
		
		for(int i=0; i < x.length; i++)
		f[i] = evalDensity(x[i],data.length);
	}
	
	/*public void evalDensity(double[] x, double[] f) {
		
		for(int i=0; i < x.length; i++)
			f[i] = evalDensity(x[i],x.length);
	}*/
	
	
}
