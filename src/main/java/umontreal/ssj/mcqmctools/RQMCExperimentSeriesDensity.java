/*
 * Class:        RQMCPointSetSeries
 * Description:  randomized quasi-Monte Carlo simulations
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
package umontreal.ssj.mcqmctools;

import umontreal.ssj.functionfit.LeastSquares;
import umontreal.ssj.hups.*;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stat.density.DEHistogram;
import umontreal.ssj.stat.density.DensityEstimator;
import umontreal.ssj.util.Chrono;
import umontreal.ssj.util.Num;
import umontreal.ssj.util.PrintfFormat;
import java.util.ArrayList;
import java.util.Iterator;

import cern.colt.Arrays;


/**
 * This class offers facilities to perform experiments to study the convergence
 * of the variance when estimating a mean (expectation) with a series of RQMC 
 * point sets usually of the same type, but different size @f$n@f.
 * The series of RQMC point sets of different sizes can be passed in an array 
 * to the constructor.   The method @f$testVarianceRate@f$ performs an experiment 
 * with a given model and the series of point sets.  One can recover the average,
 * variance, their logs in base 2, etc., in arrays, as well as the estimated 
 * linear regression of log(variance) as a function of log(n). 
 */

public class RQMCExperimentSeriesDensity extends RQMCExperimentSeries {


   /**
    * Constructor with a give series of RQMC point sets.
    *  @param theSets      the RQMC point sets
    */
   public RQMCExperimentSeriesDensity (RQMCPointSet[] theSets) {
	   super(theSets);
   }
   public RQMCExperimentSeriesDensity (ArrayList<RQMCPointSet[]> theSets) {
	   super(theSets);
   }


   /**
    * Performs an RQMC experiment with the given model, with this series of RQMC point sets and a series of density estimator.  
    * For each set in the series, computes the average, the variance, its log in base 2.
    */
   
   
   public void testVarianceRate (MonteCarloModelBounded model, int m,
			ArrayList<DensityEstimator> listDE, int numEvalPoints, 
           double[] integVariance, RQMCPointSet [] theSets) {
	   
		   for(int i=0; i<listDE.size(); i++) 		   
		   testVarianceRate (model, m, listDE.get(i), numEvalPoints, integVariance, theSets);		   
	          	  
   }
   
   /**
    * Performs an RQMC experiment with the given model, with this series of RQMC point sets.  
    * For each set in the series, computes the average, the variance, its log in base 2.
    */
   
   
   public void testVarianceRate (MonteCarloModelBounded model, int m,
			DensityEstimator DE, int numEvalPoints, 
           double[] integVariance, RQMCPointSet [] theSets) {
	int n;
	Tally statReps = new Tally();
	Chrono timer = new Chrono();
	numReplicates = m;
	this.model = model;
   if (displayExec) {
   	System.out.println("\n ============================================= ");
   	System.out.println("RQMC simulation for density estimation, for unknown density:  ");
   	System.out.println("Model: " + model.toString());
   	System.out.println(" Number of indep copies m  = " + m);
   	System.out.println(" Point sets: " + theSets[0].toString() + "\n");
	System.out.println("    n     CPU time         mean       log2(var) ");	    	
   }

    int r=8, l=2;
	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		size[s] = n;
		double[][] data = new double[m][];
		log2n[s] = Num.log2(n);	
		if ( DE == new DEHistogram(DE.getA(),DE.getB()) )
		    log2h[s] = Math.log((DE.getB()-DE.getA())/Math.pow(8, l));
			 //log2h[s] = DE.geth();
		//log2h[s] = DE.geth();		
		//log2h[s]= r*Math.log(2);
		log2h[s]= Math.log(1/Math.pow(2,r));
		l++;
		r++;
		RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statReps, data);		
		integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model, n, m, data, DE, numEvalPoints);
		mean[s] = statReps.average();
	    log2Var[s] = Num.log2(integVariance[s]);
	    if (displayExec) {
		   System.out.println("  " + n + "     " + timer.format() + 
		              "   " + PrintfFormat.f(10, 5, mean[s]) + 
				      "   " + PrintfFormat.f(7, 2, log2Var[s]));
	    }
	}	
	 
   cpuTime = timer.format();	   
}
   
   



   /**
    * Returns the vector of log_2(n).
    */
   public double[] getLog2n() {
      return log2n;
   }

   
   /**
    * Produces and returns a report on the last experiment.
    * @param numSkip  The first numSkip values of n are skipped for the regression
    * @param details  If true, gives values (mean, log variance,...) for each n.
    * @return  Report as a string.
    */
	public String reportFixedh (boolean details) {
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for density estimation, with unknown density: \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + numReplicates + "\n");
		sb.append(" Point sets: " + this.toString() + "\n\n");
		sb.append("RQMC integrated variance (IV) \n");
		if (details) {
			sb.append("    n      mean       log2(var) \n");
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				sb.append(" " + size[s] + " " + PrintfFormat.f(10, 5, mean[s]) +
				          " " + PrintfFormat.f(7, 2, log2Var[s]) + "\n");
			}
		}
		double[] regCoeff = regressionLogVariance (numSkipRegression);
		sb.append("  Slope of log2(var) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		sb.append("    constant term      = " + PrintfFormat.f(8, 5, Math.exp(regCoeff[0])) + "\n\n");
		sb.append("  Total CPU Time = " + cpuTime + "\n");
		sb.append("-----------------------------------------------------\n");		
		return sb.toString();
	}
	
	
	public double[] regressionLogVarianceVariedh (int numSkip) {
		double[][] x2 = new double[numSets-numSkip][2];
		double [] y2 = new double[numSets-numSkip];
		for (int i = 0; i < numSets-numSkip; ++i) {
			x2[i][0] = log2n[i+numSkip];
			x2[i][1] = log2h[i+numSkip];			
			y2[i] = log2Var[i+numSkip];
		}
		return LeastSquares.calcCoefficients0(x2, y2);
	}
	
	public String reportVariedH (boolean details, double alpha) {
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for density estimation, with unknown density: \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + numReplicates + "\n");
		sb.append(" Point sets: " + this.toString() + "\n\n");
		sb.append("RQMC integrated variance (IV) \n");
		if (details) {
			sb.append("    n      mean       log2(var) \n");
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				sb.append(" " + size[s] + " " + PrintfFormat.f(10, 5, mean[s]) +
				          " " + PrintfFormat.f(7, 2, log2Var[s]) + "\n");
			}
		}
		double[] regCoeff = regressionLogVarianceVariedh (numSkipRegression);
		//sb.append("  Slope of log2(var) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		//sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");
		sb.append("  C     = " + Math.exp(regCoeff[0]) + "\n");
		sb.append("  Slope of log2(var) : beta  = " + -regCoeff[1] + "\n");
		sb.append("  delta = " + -regCoeff[2] + "\n");		
		//sb.append("  gamma = " + (-regCoeff[1])/(alpha - regCoeff[2]) + "\n");	
		//sb.append("  nu    = " + (-alpha * regCoeff[1])/(alpha - regCoeff[2]) + "\n\n");	
		sb.append("  Total CPU Time = " + cpuTime + "\n");
		sb.append("-----------------------------------------------------\n");		
		return sb.toString();
	}
	
	
	
	/**
	 * Performs an experiment (testVarianceRate)  with the given model for each point set series in the given list,
	 * and returns a report as a string. 
	 * 
	 */
	
	public String TestRQMCManyPointTypes (MonteCarloModelBounded model, 
			ArrayList<RQMCPointSet[]> list, int m,
			ArrayList<DensityEstimator> listDE, int numEvalPoints, 
            boolean details) {
		StringBuffer sb = new StringBuffer("");
		numReplicates = m;	
		double[] integVariance= new double[numSets];   // Will contain the IV estimates, for each n.
		for(int i=0; i < listDE.size(); i++) {
		  for (RQMCPointSet[] ptSeries : list) {			
         	testVarianceRate (model, m, listDE.get(i), numEvalPoints, integVariance, ptSeries);
			sb.append (reportFixedh (details));			
		  }
	    }
		return sb.toString();
	}
	
	
	/**
	 * Performs an experiment (testVarianceRate) for each point set series in the given list,
	 * and returns a report as a string. 
	 * 
	 */
	
	public String toString () {
		return theSets[0].toString();
	}
}