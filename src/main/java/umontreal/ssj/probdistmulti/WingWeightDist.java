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
public class WingWeightDist extends ContinuousDistributionMulti {
   protected int dim;
   protected double[] x;
 
   
	   public WingWeightDist (double [] x ) {
		   this.x=x;
		   this.dim = x.length;
   }



	   
	 public double density (double[] x) {   
		 double Sw = x[0], Wfw  = x[1], A = x[2], LamCaps = x[3], q = x[4], lam = x[5], tc = x[6], Nz = x[7], Wdg = x[8], Wp =x[9];	   		 

		double fact1 = 0.036 * Math.pow(Sw, 0.758) * Math.pow(Wfw, 0.0035);
		double fact2 = Math.pow(A / (Math.pow(Math.cos(LamCaps),2)),0.6);
		double fact3 = Math.pow(q,0.006) * Math.pow(lam, 0.04);
		double fact4 = Math.pow(100*tc / Math.cos(LamCaps),(-0.3));
		double fact5 = Math.pow(Nz*Wdg, 0.49);

		double term1 = Sw * Wp;

		 return fact1*fact2*fact3*fact4*fact5 + term1;
	       
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