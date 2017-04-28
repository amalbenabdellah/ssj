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
 * its value at @f$m@f$ selected points. For details on how the kernel
 * density is defined, and how to select the kernel and the bandwidth
 * @f$h@f$, see the documentation of class
 * @ref umontreal.ssj.randvar.KernelDensityGen in package `randvar`.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class DEKernelDensity implements DensityEstimator {
	double a, b, h; 
	KernelDensity kerD;
	EmpiricalDist dist;
	ContinuousDistribution kern ;
	
	public DEKernelDensity(double a, double b, double h, ContinuousDistribution kern  ) {
		
		this.a = a;
		this.b = b;
		this.h = h;
		//this.dist = dist;
		this.kern = kern;
	}
	
    public void setRange(double a, double b) {
		
		this.a = a;
		this.b = b;
	}
	
	public void constructDensity(double[] x) {
		
		kerD = new KernelDensity();
		dist= new EmpiricalDist(x);
		
		
	}

	public double evalDensity(double x, int n) {
		return kerD.estimate(dist, kern, h, x);
		
	}

	public void evalDensity(double[] x, double[] f, double[] data ) {
		
		for(int i=0; i < x.length; i++)
			f[i] = evalDensity(x[i], data.length);
		
		//f = kerD.computeDensity(dist, kern, h, x);
		
	}
	

   /**
    * Given the empirical distribution `dist`, this method computes the
    * kernel density estimate at each of the @f$m@f$ points
    * <tt>Y[</tt>@f$j@f$<tt>]</tt>, @f$j= 0, 1, …, (m-1)@f$, where @f$m@f$
    * is the length of `Y`, the kernel is `kern.density(x)`, and the
    * bandwidth is @f$h@f$. Returns the estimates as an array of @f$m@f$
    * values.
    */
  /* public static double[] computeDensity (double[] data,
                                          ContinuousDistribution kern,
                                          double h, double[] Y) {
      int m = Y.length;
      double[] u = new double[m];
      for (int j = 0; j < m; j++)
         u[j] = estimate (dist, kern, h, Y[j]);
      return u;
   }*/

   /**
    * Similar to method
    * #computeDensity(EmpiricalDist,ContinuousDistribution,double,double[])
    * above
    * , but the bandwidth @f$h@f$ is obtained from the method
    * {@link umontreal.ssj.randvar.KernelDensityGen.getBaseBandwidth(EmpiricalDist)
    * getBaseBandwidth(dist)} in package `randvar`.
    */
   /*public static double[] computeDensity (EmpiricalDist dist,
                                          ContinuousDistribution kern,
                                          double[] Y) {
      double h = KernelDensityGen.getBaseBandwidth(dist);
      return computeDensity (dist, kern, h, Y);
   }*/

public double geth() {
	// TODO Auto-generated method stub
	return h;
}



public double getA() {
	// TODO Auto-generated method stub
	return a;
}

public double getB() {
	// TODO Auto-generated method stub
	return b;
}

}