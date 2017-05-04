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
public class PistonDist extends ContinuousDistributionMulti {
   protected int dim;
   protected double[] x;
 
   
	   public PistonDist (double [] x ) {
		   this.x=x;
		   this.dim = x.length;
   }



	   
	 public double density (double[] x) {      
	      double M = x[0], S  = x[1], V0 = x[2], k = x[3], P0 = x[4], Ta = x[5], T0 = x[6];	     	 

	    double Aterm1 = P0 * S;
	    double Aterm2 = 19.62 * M;
	    double Aterm3 = -k*V0 / S;
	    double A = Aterm1 + Aterm2 + Aterm3;

	    double Vfact1 = S / (2*k);
	    double Vfact2 = Math.sqrt(Math.pow(A, 2) + 4*k*(P0*V0/T0)*Ta);
	    double V = Vfact1 * (Vfact2 - A);

	    double fact1 = M;
	    double fact2 = k + (Math.pow(S, 2))*(P0*V0/T0)*(Ta/(Math.pow(V, 2)));

	     return 2 * Math.PI * Math.sqrt(fact1/fact2);


	      
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