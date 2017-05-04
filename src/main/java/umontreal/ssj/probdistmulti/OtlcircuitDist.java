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
public class OtlcircuitDist extends ContinuousDistributionMulti {
   protected int dim;
   protected double[] x;
 
   
	   public OtlcircuitDist (double [] x ) {
		   this.x=x;
		   this.dim = x.length;
   }



	   
	 public double density (double[] x) {      
	      double Rb1 = x[0], Rb2  = x[1], Rf = x[2], Rc1 = x[3], Rc2 = x[4], beta = x[5];
	      

	     double  Vb1 = 12*Rb2 / (Rb1+Rb2);
	     double  term1a = (Vb1+0.74) * beta * (Rc2+9);
	     double  term1b = beta*(Rc2+9) + Rf;
	     double  term1 = term1a / term1b;

	     double term2a = 11.35 * Rf;
	     double term2b = beta*(Rc2+9) + Rf;
	     double term2 = term2a / term2b;

	     double term3a = 0.74 * Rf * beta * (Rc2+9);
	     double term3b = (beta*(Rc2+9)+Rf) * Rc1;
	     double term3 = term3a / term3b;

	      return term1 + term2 + term3;
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