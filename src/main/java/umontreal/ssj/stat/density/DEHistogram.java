package umontreal.ssj.stat.density;

import umontreal.ssj.stat.HistogramOnly;
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
	//TallyHistogram hist;
	ScaledHistogram histDensity;

	DEHistogram (double a, double b) {
		this.a=a;
		this.b=b;
		
	}
    DEHistogram (double a, double b,int  m) {
    	this.a=a;
		this.b=b;
		this.m=m;
	}

	
	public void setRange(double a, double b) {
		// TODO Auto-generated method stub
		this.a=a;
		this.b=b;
		histDensity = null;
	}

	
	public void constructDensity(double[] x) {
		// TODO Auto-generated method stub
//		double h=(b-a)/m;
//		x[0]=a-0.5*h;
//		for(int j=1; j < m; j++)		
//			x[j]=x[j-1]+h;
		TallyHistogram hist= new TallyHistogram(a,b,m);
		hist.fillFromArray(x);

		histDensity=new ScaledHistogram(hist,1.0);
		
		//x=histDensity.getHeights();
		//for(int j=0; j < x.length; j++)
		//	x[j]=histDensity.getHeights()[j];
		
	}

	
	public double evalDensity(double x) {
		// TODO Auto-generated method stub
		//double h=(b-a)/m;
		//double[] count=histDensity.getHeights();
		
		//double evaldens=histDensity.getHeights()[(int) ((x -a) /h)];
		double evaldens=histDensity.getHeight(x);
		return evaldens;
	}

	
	public void evalDensity(double[] x, double[] f) {
		// TODO Auto-generated method stub
		for(int i=0; i < x.length; i++)
			f[i]=evalDensity(x[i]);
	}
	
}
