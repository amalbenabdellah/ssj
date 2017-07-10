/*
 * Class:        KernelDensity
 * Description:  Kernel density estimators
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since
 *
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
package umontreal.ssj.stat.density;
   import umontreal.ssj.gof.KernelDensity;
import umontreal.ssj.probdist.*;
import umontreal.ssj.randvar.KernelDensityGen;
import umontreal.ssj.stat.KernelDensityEstimator1d;

/**
 * This class provides methods to compute a kernel density estimator from a
 * set of @f$n@f$ individual observations @f$x_0, …, x_{n-1}@f$, and returns
 * its value at a  set of selected points. 
 */
public class DEKernelDensity implements DensityEstimator {
	double a, b, h; 
	KernelDensity kerD;
	EmpiricalDist dist;
	ContinuousDistribution kern ;
	
public DEKernelDensity(double a, double b, double h ) {
		
		this.a = a;
		this.b = b;
		this.h = h;
		
	}
	
	public DEKernelDensity(double a, double b, double h, ContinuousDistribution kern  ) {
		
		this.a = a;
		this.b = b;
		this.h = h;
		this.kern = kern;
	}
public DEKernelDensity(double a, double b ) {
		
		this.a = a;
		this.b = b;
	}
	
    public void setRange(double a, double b) {
		
		this.a = a;
		this.b = b;
	}
	
    
    /**	
     * Constructs a density estimator from 
     * the data points in vector x.
     * @param x  data points.
     */
	public void constructDensity(double[] x) {
		
		kerD = new KernelDensity();
		dist= new EmpiricalDist(x);
		
		
	}

	/**
	 * Returns the value of the density evaluated at x.
	 */
	public double evalDensity(double x, int n) {
		return kerD.estimate(dist, kern, h, x);
		
	}

	/**	
     * Returns in array density the value of the density at the evaluation points in x.
     * These two arrays must have the same size. 
     * @param x  evaluation points
     * @param density  values of the density
     * @param data array of observations
     */
	public void evalDensity(double[] x, double[] f, double[] data ) {
		
		for(int i=0; i < x.length; i++)
			f[i] = evalDensity(x[i], data.length);
		
	}
	
	public String toString(){
		return "KDE Gaussian_"+ geth();
	}
	

   

	
	public void seth (double h){
		this.h=h;
	}
	
	
public double geth() {
	
	return h;
}



public double getA() {
	
	return a;
}

public double getB() {
	
	return b;
}
public void setM(int m){
		
};
public int getM(){
	return (int)(b-a/h);
};

}