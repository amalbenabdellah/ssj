/*
 * Class:        MultiNormalDist
 * Description:  multinormal distribution
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
package umontreal.ssj.probdistmulti;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;

/**
 * Implements the abstract class  @ref ContinuousDistributionMulti for the
 * *multinormal* distribution with mean vector @f$\boldsymbol{\mu}@f$ and
 * covariance matrix @f$\boldsymbol{\Sigma}@f$. The probability density is
 * @anchor REF_probdistmulti_MultiNormalDist_eq_fMultinormal
 * @f[
 *   f(\mathbf{x}) = \frac{1}{\sqrt{(2\pi)^d \det\boldsymbol{\Sigma}}} \exp\left(-\frac{1}{2}(\mathbf{x}- \boldsymbol{\mu})^T \boldsymbol{\Sigma}^{-1} (\mathbf{x}- \boldsymbol{\mu})\right) \tag{fMultinormal}
 * @f]
 * where @f$\mathbf{x}= (x_1,…,x_d)@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdistmulti_continuous
 */
public class RobotDist extends ContinuousDistributionMulti {
   protected int dim;
   protected double[] x;
 
   
	   public RobotDist (double [] x ) {
		   this.x=x;
		   this.dim = x.length;
   }



	   
	 public double density (double[] x) {      
	       double [] theta = new double [4], L= new double [4], Li= new double [4], thetai= new double [4];
	      for ( int i = 0; i < 4; i++)
	    	  theta[i] = x[i];
	      for ( int i = 4; i < 8; i++)
	    	   L[i-4] = x[i];

	   double   sumu = 0, sumv = 0;
	   for ( int ii = 0; ii < 4; ii++){
	          Li[ii] = L[ii];
	        double  sumtheta = 0;
	      for ( int jj = 0; jj <ii; jj++){	          
	              thetai[jj] = theta[jj];
	              sumtheta = sumtheta + thetai[jj];
	      }
	          sumu = sumu + Li[ii] * Math.cos(sumtheta);
	          sumv = sumv + Li[ii]*Math.sin(sumtheta);
	   }

	      return Math.pow((Math.pow(sumu, 2) + Math.pow(sumv, 2)),(0.5));
	 }

	 
	 public int getDimension() {
	      return dim;
	   }

	public double[] getMean() {
		// TODO Auto-generated method stub
		return null;
	}



	public double[][] getCovariance() {
		// TODO Auto-generated method stub
		return null;
	}



	public double[][] getCorrelation() {
		// TODO Auto-generated method stub
		return null;
	}

   




}