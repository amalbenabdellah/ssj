package umontreal.ssj.stat.density;

import umontreal.ssj.stat.AveragedShiftedHistogram;
import umontreal.ssj.stat.ScaledHistogram;
import umontreal.ssj.stat.TallyHistogram;

/**
 *
* An averaged shifted histogram used a a density estimator.
* 
*/
public class  DEAveragedShiftedHistogram  implements DensityEstimator {
	double a,b;
	int m, r;
	ScaledHistogram histDensity;
	ScaledHistogram histAsh;
	AveragedShiftedHistogram ASH;
	AveragedShiftedHistogram ASH2;
	
	
	DEAveragedShiftedHistogram (double a, double b) {
		this.a = a;
		this.b = b;
		histDensity = null;
		histAsh = null;
		ASH = null;
		ASH2 = null;
		
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
	public void seth(double h) {
		this.m = (int) ((b-a)/h);
	}
	public double geth() {
		return (b-a)/m;
	}
	public void setM(int m) {
		this.m = m;
	}
	public int getR() {
		return r;
	}
	public void setR(int r) {
		this.r = r;
	}
	public ScaledHistogram getHistDensity() {
		return histDensity;
	}
	public void setHistDensity(ScaledHistogram histDensity) {
		this.histDensity = histDensity;
	}
	public ScaledHistogram getHistAsh() {
		return histAsh;
	}
	public void setHistAsh(ScaledHistogram histAsh) {
		this.histAsh = histAsh;
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

	
	// Estimate the density 
	
	public void constructDensity(double[] x) {
		
		/*TallyHistogram hist = new TallyHistogram(a,b,m);
		hist.fillFromArray(x);
		
		histDensity=new ScaledHistogram(hist,1.0);
		ScaledHistogram histAsh=histDensity.averageShiftedHistogramTrunc(r);
		histAsh.rescale(1.0);*/
		
		
		TallyHistogram hist = new TallyHistogram(a,b,m*r);
	    hist.fillFromArray(x);
		ASH = new AveragedShiftedHistogram(hist,1.0);
		ASH2 = ASH.averageShiftedHistogramTrunc(r);
		ASH2.rescale(1.0);
	}

	// Evaluate  the density  at a point
	
	public double evalDensity(double x, int n) {
		double h = (b - a)/ m;	
		double evaldens = ASH2.getHeights()[(int) ((x - a)*r /h)]/(n*h*r*r);
		//double evaldens = histAsh.getHeight(x);
		return evaldens;
	}
	
  
	
	/*public void evalDensity(double[] x, double[] f) {
		
		for(int i=0; i < x.length; i++)
			f[i] = evalDensity(x[i], x.length);
	}*/
	
	  // Evaluate the density at eval points
	
    public void evalDensity(double[] x, double[] f, double[] data) {
    	
    	for(int i=0; i < x.length; i++)
			f[i] = evalDensity(x[i], data.length);
	}
	
}
