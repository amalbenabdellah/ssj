package umontreal.ssj.stat.density;

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
	
	
	DEAveragedShiftedHistogram (double a, double b) {
		this.a=a;
		this.b=b;
		histDensity=null;
		histAsh=null;
		
	}
	DEAveragedShiftedHistogram (double a, double b, int m, int r) {
	this.a=a;
	this.b=b;
	this.m=m;
	this.r=r;
	histDensity=null;
	histAsh=null;

	}
	public void setRange(double a, double b) {
		// TODO Auto-generated method stub
		this.a=a;
		this.b=b;
		histDensity=null;
		histAsh=null;
	}

	
	public void constructDensity(double[] x) {
		// TODO Auto-generated method stub
		TallyHistogram hist= new TallyHistogram(a,b,m);
		hist.fillFromArray(x);
		histDensity=new ScaledHistogram(hist,1.0);
		ScaledHistogram histAsh=histDensity.averageShiftedHistogramTrunc(r);
		histAsh.rescale(1.0);
	}

	
	public double evalDensity(double x) {
		// TODO Auto-generated method stub
		double evaldens=histAsh.getHeight(x);
		return evaldens;
	}

	
	public void evalDensity(double[] x, double[] f) {
		// TODO Auto-generated method stub
		for(int i=0; i < x.length; i++)
			f[i]=evalDensity(x[i]);
	}
	
}
