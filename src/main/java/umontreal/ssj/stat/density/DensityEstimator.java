/**
 * 
 */
package umontreal.ssj.stat.density;

/**
 * @author Lecuyer
 * 
 * An interface for a density estimator. 
 *
 */
public interface DensityEstimator {
	
	
	/**
	 * From now on, the density will be estimated over the interval [a,b].
	 * 
	 */
	public void setRange (double a, double b);
	 
	 
    /**	
     * Constructs a density estimator from 
     * the data points in vector x.
     * @param x  data points.
     */
	public void constructDensity (double[] x);
	

	/**
	 * Returns the value of the density evaluated at x.
	 */
	public double evalDensity (double x, int n);

	
    /**	
     * Returns in array density the value of the density at the evaluation points in x.
     * These two arrays must have the same size. 
     * @param x  evaluation points
     * @param density  values of the density
     */
	
	public void evalDensity (double[] x, double[] density, double[] data);
	/**	
	 *  if the density is kernel h is the bandwidth else if the density is a histogram  h represents (rang/numbins)
	 */
	
	public double geth();
	public void seth (double h);
	public double getA();
	public double getB();
	public void setM(int m);
	public int getM();

		
	}
		


