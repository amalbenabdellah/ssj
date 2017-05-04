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
public class BoreholeDist extends ContinuousDistributionMulti {
   protected int dim;
   protected double[] x;
   
   /*public BoreholeDist (double rw , double r, double Tu, double Hu, double Tl, double Hl, double L, double Kw ) {
	   
   }*/
	   public BoreholeDist (double [] x ) {
		   this.x=x;
   }

	   public double density (  double rw , double r, double Tu, double Hu, double Tl, double Hl, double L, double Kw){

      double frac1 = 2 * Math.PI * Tu * (Hu-Hl);

      double frac2a = 2*L*Tu / (Math.log(r/rw)*Math.pow(rw, 2)*Kw);
      double frac2b = Tu / Tl;
      double frac2 = Math.log(r/rw) * (1+frac2a+frac2b);

      return frac1 / frac2;
            
   }
	   
	 public double density (double[] x) {      
	      double rw = x[0], r  = x[1], Tu = x[2], Hu = x[3], Tl = x[4], Hl = x[5], L  = x[6], Kw = x[7];
	      double frac1 = 2 * Math.PI * Tu * (Hu-Hl);

	      double frac2a = 2*L*Tu / (Math.log(r/rw)*Math.pow(rw, 2)*Kw);
	      double frac2b = Tu / Tl;
	      double frac2 = Math.log(r/rw) * (1+frac2a+frac2b);

	      return frac1 / frac2;
	 }


	 public int getDimension() {
	      return dim;
	   }
	@Override
	public double[] getMean() {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public double[][] getCovariance() {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public double[][] getCorrelation() {
		// TODO Auto-generated method stub
		return null;
	}

   




}