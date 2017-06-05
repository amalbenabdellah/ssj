package umontreal.ssj.mcqmctools;
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


import umontreal.ssj.charts.XYLineChart;
import umontreal.ssj.functionfit.LeastSquares;
import umontreal.ssj.hups.*;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stat.density.DEAveragedShiftedHistogram;
import umontreal.ssj.stat.density.DEHistogram;
import umontreal.ssj.stat.density.DensityEstimator;
import umontreal.ssj.util.Chrono;
import umontreal.ssj.util.Num;
import umontreal.ssj.util.PrintfFormat;

import java.io.FileWriter;
import java.io.IOException;
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


	XYLineChart chartIV = new XYLineChart ();   
	double[] log2h; // log_2 of h
	double[] log2IV;  // log_2 of Integrated Variance
	
	/**
	    * Constructor with a give series of RQMC point sets.
	    *  @param theSets      the RQMC point sets
	    */
	
   public RQMCExperimentSeriesDensity (RQMCPointSet[] theSets) {
	   super(theSets);
	  
	   log2h = new double[numSets];  	   
	   log2IV = new double[numSets];
	   	   
   }
 
   /**
    * Performs an RQMC experiment with the given model, with this series of RQMC point sets for fixed h.  
    * For each set in the series, computes the average, the variance, the integrated variance, its log in base 2.
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
	System.out.println("    n     CPU time    log2Var   log2(IV) ");	    	
   }


	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		size[s] = n;
		double[][] data = new double[m][];
		log2n[s] = Num.log2(n);	
		RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statReps, data);	
		mean[s] = statReps.average();
		variance [s] = statReps.variance();
		log2Var[s] = Num.log2(variance[s]);	
		integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model,  m, data, DE, numEvalPoints);
	    log2IV[s] = Num.log2(integVariance[s]);
	    if (displayExec) {
		   System.out.println("  " + n + "     " + timer.format()+ 
				      "   " + PrintfFormat.f(7, 2, log2Var[s]) + 
				      "   " + PrintfFormat.f(7, 2, log2IV[s]));
	    }
	}	
	 
   cpuTime = timer.format();	   
}
   
   /**
    * Performs an RQMC experiment with the given model, with this series of RQMC point sets and with varied h and varied n .  
    * For each set in the series, computes the average, the variance, the integrated variance, its log in base 2.
    */
   
   
   public void testVarianceRateVariedhn (MonteCarloModelBounded model, int m,
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
	System.out.println("    n     CPU time     log2Var log2(IV) ");	    	
   }

    double  l=1;
	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		size[s] = n;
		double[][] data = new double[m][];
		
		log2n[s] = Num.log2(n);	
		
	
		
		if ( DE == new DEHistogram(DE.getA(),DE.getB()) ){
		    log2h[s] = Math.log((DE.getB()-DE.getA())/Math.pow(4, l));
		    DE.seth((DE.getB()-DE.getA())/Math.pow(4, l));
		    l++;
			}
		else if ( DE == new DEAveragedShiftedHistogram(DE.getA(),DE.getB()) ){
			log2h[s] = Math.log((DE.getB()-DE.getA())/(32*Math.pow(4, l)));
		    DE.seth((DE.getB()-DE.getA())/(32*Math.pow(4, l)));	
		    l++;
			
		}
		else {
			log2h[s] =- Math.log(Math.pow(4, l));
		}

		RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statReps, data);		
		integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model,  m, data, DE, numEvalPoints);
		mean[s] = statReps.average();
		variance [s] = statReps.variance();
		log2Var[s] = Num.log2(variance[s]);	
	    log2IV[s] = Num.log2(integVariance[s]);
	    if (displayExec) {
	    	System.out.println("  " + n + "     " + timer.format()+ 
				      "   " + PrintfFormat.f(7, 2, log2Var[s]) + 
				      "   " + PrintfFormat.f(7, 2, log2IV[s]));
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


	
   /** Fits a regression starting from numSkip,  between the arrays x, y and z and return the coefficients
    */
   public double[] slope ( double[] x, double [] y,  double []z, int numSkip) {
		double[][] x2 = new double[numSets-numSkip][2];
		double [] y2 = new double[numSets-numSkip];
		for (int i = 0; i < numSets-numSkip; ++i) {
			x2[i][0] = x[i+numSkip];
			x2[i][1] = y[i+numSkip];			
			y2[i] = z[i+numSkip];
		}
		return LeastSquares.calcCoefficients0(x2, y2);
	}
   
   
   /** Fits a regression starting from numSkip,  between the arrays x and y and return the coefficients
    */
  
  public double[] slope (double[] x, double [] y, int numSkip) {
		double[] x2 = new double[numSets-numSkip];
		double [] y2 = new double[numSets-numSkip];
		for (int i = 0; i < numSets-numSkip; ++i) {
			x2[i] = x[i+numSkip];	
			y2[i] = y[i+numSkip];
		}
		return LeastSquares.calcCoefficients(x2, y2,1);
	}
	
	
	
   
   /**
    * Produces and returns a report on the last experiment.
    * @param numSkip  The first numSkip values of n are skipped for the regression
    * @param details  If true, gives values (mean, log variance,...) for each n.
    * @return  Report as a string.
    */
	
	
	
	public String reportVarianceVariedhn (boolean details, String densityestimator) {
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for density estimation, with unknown density: \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + numReplicates + "\n");
		sb.append(" Density Esimator  = " + densityestimator + "\n");
		sb.append(" Point sets: " + this.toString() + "\n\n");
		sb.append("RQMC integrated variance (IV) \n");
		if (details) {
			sb.append("    n    log2Var log2(IV) \n");
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				sb.append(" " + size[s] +  " " + PrintfFormat.f(7, 2, log2Var[s]) +
				          " " + PrintfFormat.f(7, 2, log2IV[s]) + "\n");
			}
		}
		double[] regCoeff = slope (log2n, log2h, log2IV, numSkipRegression);
		sb.append("  C     = " + Math.exp(regCoeff[0]) + "\n");
		sb.append("  Slope of log2(IV) : beta  = " + -regCoeff[1] + "\n");
		sb.append("  delta = " + -regCoeff[2] + "\n");		
		regCoeff = slope(log2n, log2Var, numSkipRegression);
		sb.append("  Slope of log2(Var) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		sb.append("  Total CPU Time = " + cpuTime + "\n");
		sb.append("-----------------------------------------------------\n");		
		return sb.toString();
	}
	
	public String reportD(boolean details , String densityestimator) {
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for density estimation: \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + numReplicates + "\n");
		sb.append(" Density Esimator  = " + densityestimator + "\n");
		sb.append(" Point sets: " + this.toString() + "\n\n");
		sb.append("RQMC integrated Variance IV \n");
		if (details) {
			sb.append("    n    log2Var log2(IV) \n");
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				sb.append(" " + size[s] +  " " + PrintfFormat.f(7, 2, log2Var[s]) +
				          " " + PrintfFormat.f(7, 2, log2IV[s]) + "\n");
			}
		}
		double[] regCoeff = slope (log2n, log2IV, numSkipRegression);
		
		sb.append("  Slope of log2(IV) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");		
		regCoeff = slope(log2n, log2Var, numSkipRegression);
		sb.append("  Slope of log2(Var) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		sb.append("  Total CPU Time = " + cpuTime + "\n");
		sb.append("-----------------------------------------------------\n");		
		return sb.toString();
	}
	
	
	/**
	 * Performs an experiment (testIntegratedVarianceRate) for fixed h  with the given model for each point set series in the given list,
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
         	this.testVarianceRate (model, m, listDE.get(i), numEvalPoints, integVariance, ptSeries);
			sb.append (reportD (details, listDE.get(i).toString()));	
			makePlotsVar (numSets, m, (model.toString()).split(" ")[0], listDE.get(i).toString(), ptSeries.toString());
		  }
	    }
		return sb.toString();
	}
	
	
	/**
	 * Performs an experiment (testIntegratedVarianceRate) for varied  h and n  with the given model for each point set series in the given list,
	 * and returns a report as a string. 
	 * 
	 */
	public String TestRQMCManyPointTypesVariedhn (MonteCarloModelBounded model, 
			ArrayList<RQMCPointSet[]> list, int m,
			ArrayList<DensityEstimator> listDE, int numEvalPoints, 
            boolean details) {
		StringBuffer sb = new StringBuffer("");
		numReplicates = m;	
		double[] integVariance= new double[numSets];   // Will contain the IV estimates, for each n.
		for(int i=0; i < listDE.size(); i++) {
		  for (RQMCPointSet[] ptSeries : list) {	
         	this.testVarianceRateVariedhn (model, m, listDE.get(i), numEvalPoints, integVariance, ptSeries);
			sb.append (reportVarianceVariedhn (details, listDE.get(i).toString()));	
			makePlotsVar (numSets, m, (model.toString()).split(" ")[0], listDE.get(i).toString(), ptSeries.toString());
		  }
	    }
		return sb.toString();
	}
	
	
	/** Method that produces Latex plots for the experiment for the IV
	 */
	public void makePlotsVar (int numSets, int m, String descModel, String densityEstimator, String descPoints) {
		
		try {
			    chartIV.add(log2n, log2IV, densityEstimator, " ");
			FileWriter file = new FileWriter(descModel + "_" +  densityEstimator + "_" + descPoints + "_IV.tex");
			file.write(chartIV.toLatex(12, 8));
			file.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	

	
	public String toString () {
		return theSets[0].toString();
	}
}