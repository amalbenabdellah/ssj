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


	XYLineChart chart = new XYLineChart ();
   /**
    * Constructor with a give series of RQMC point sets.
    *  @param theSets      the RQMC point sets
    */
	double[] log2h; // log_2 of h
	public double[][] meanD ;
	public double[][] log2VarS; 
	double[] log2IV;
	
	
   public RQMCExperimentSeriesDensity (RQMCPointSet[] theSets) {
	   super(theSets);
	  
	   log2h = new double[numSets];  
	   
	   log2IV = new double[numSets];
	   
	   
   }
   public RQMCExperimentSeriesDensity (ArrayList<RQMCPointSet[]> theSets) {
	   super(theSets);
	   log2h = new double[numSets];  
	  
	   log2IV = new double[numSets];
   }


   /**
    * Performs an RQMC experiment with the given model, with this series of RQMC point sets and a series of density estimator.  
    * For each set in the series, computes the average, the variance, its log in base 2.
    */
   
   public void testVarianceRate (MonteCarloModelBounded model, int m,
		   ArrayList<DensityEstimator> listDE, int numEvalPoints, 
           double[][] integVariance, RQMCPointSet [] theSets) {
	int n;
	Tally statReps = new Tally();
	Tally statKS = new Tally();
	Tally statCVM = new Tally();
	Chrono timer = new Chrono();
	numReplicates = m;
	this.model = model;
   if (displayExec) {
   	System.out.println("\n ============================================= ");
   	System.out.println("RQMC simulation for density estimation, for unknown density:  ");
   	System.out.println("Model: " + model.toString());
   	System.out.println(" Number of indep copies m  = " + m);
   	System.out.println(" Point sets: " + theSets[0].toString() + "\n");
	System.out.println("    n     CPU time         mean       log2(IV) ");	    	
   }


    
	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		for(int i=0; i<listDE.size(); i++){	
		size[s] = n;
		double[][] data = new double[m][];
		
		log2n[s] = Num.log2(n);	
		meanD= new double[listDE.size()][numSets];
		log2VarS= new double [listDE.size()][numSets];
		RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statReps, data);		
		integVariance[i][s]=RQMCExperimentDensity.computeDensityVariance (model, m, data, listDE.get(i), numEvalPoints);
		meanD[i][s] = statReps.average();
	    log2VarS[i][s] = Num.log2(integVariance[i][s]);
	    if (displayExec) {
	    	
	    	if (displayExec) {
	 		   System.out.println("  " + n + "     " + timer.format() + PrintfFormat.f(7, 2,meanD[i][s])+
	 				      "   " + PrintfFormat.f(7, 2, log2VarS[i][s]));
	    	
	    }
	}	}
    }
   cpuTime = timer.format();	   
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
	Tally statKS = new Tally();
	Tally statCVM = new Tally();
	Chrono timer = new Chrono();
	numReplicates = m;
	this.model = model;
   if (displayExec) {
   	System.out.println("\n ============================================= ");
   	System.out.println("RQMC simulation for density estimation, for unknown density:  ");
   	System.out.println("Model: " + model.toString());
   	System.out.println(" Number of indep copies m  = " + m);
   	System.out.println(" Point sets: " + theSets[0].toString() + "\n");
	System.out.println("    n     CPU time    Var  KS CVM    log2(IV) ");	    	
   }


	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		size[s] = n;
		double[][] data = new double[m][];
		log2n[s] = Num.log2(n);	
		RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statReps, data);	
		mean[s] = statReps.average();
		variance [s] = statReps.variance();
		//KS[s] = statKS.average();
		//CVM[s] = statCVM.average();
		integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model,  m, data, DE, numEvalPoints);
	    log2IV[s] = Num.log2(integVariance[s]);
	    if (displayExec) {
		   System.out.println("  " + n + "     " + timer.format() + 
				      "   " + PrintfFormat.f(7, 2, log2IV[s]));
	    }
	}	
	 
   cpuTime = timer.format();	   
}
   
   
   public void testVarianceRateVariedhn (MonteCarloModelBounded model, int m,
			DensityEstimator DE, int numEvalPoints, 
           double[] integVariance, RQMCPointSet [] theSets) {
	int n;
	Tally statReps = new Tally();
	Tally statKS = new Tally();
	Tally statCVM = new Tally();
	Chrono timer = new Chrono();
	numReplicates = m;
	this.model = model;
   if (displayExec) {
   	System.out.println("\n ============================================= ");
   	System.out.println("RQMC simulation for density estimation, for unknown density:  ");
   	System.out.println("Model: " + model.toString());
   	System.out.println(" Number of indep copies m  = " + m);
   	System.out.println(" Point sets: " + theSets[0].toString() + "\n");
	System.out.println("    n     CPU time      log2(IV) ");	    	
   }

    double  l=1;
    		// t=-5;
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
			l=l+0.1;
			/*log2h[s] =t*Math.log(2);
		    DE.seth(t* Math.log(2));	
		    t++;*/
		}
			
		
		
	
		RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statReps, data);		
		integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model,  m, data, DE, numEvalPoints);
		
		mean[s] = statReps.average();
		variance [s] = statReps.variance();
		//KS[s] = statKS.average();
		//CVM[s] = statCVM.average();
	    log2IV[s] = Num.log2(integVariance[s]);
	    if (displayExec) {
		   System.out.println("  " + n + "     " + timer.format() + 
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

   
   /**
    * Produces and returns a report on the last experiment.
    * @param numSkip  The first numSkip values of n are skipped for the regression
    * @param details  If true, gives values (mean, log variance,...) for each n.
    * @return  Report as a string.
    */
	/*public String reportVarianceFixedh (boolean details) {
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
	}*/
	
	
	public double[] regressionLogVarianceVariedhn (int numSkip) {
		double[][] x2 = new double[numSets-numSkip][2];
		double [] y2 = new double[numSets-numSkip];
		for (int i = 0; i < numSets-numSkip; ++i) {
			x2[i][0] = log2n[i+numSkip];
			x2[i][1] = log2h[i+numSkip];			
			y2[i] = log2IV[i+numSkip];
		}
		return LeastSquares.calcCoefficients0(x2, y2);
	}
	
	
   
   /**
    * Produces and returns a report on the last experiment.
    * @param numSkip  The first numSkip values of n are skipped for the regression
    * @param details  If true, gives values (mean, log variance,...) for each n.
    * @return  Report as a string.
    */
	
	
	
	public String reportVarianceVariedhn (boolean details) {
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for density estimation, with unknown density: \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + numReplicates + "\n");
		sb.append(" Point sets: " + this.toString() + "\n\n");
		sb.append("RQMC integrated variance (IV) \n");
		if (details) {
			sb.append("    n     log2(IV) \n");
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				sb.append(" " + size[s] + 
				          " " + PrintfFormat.f(7, 2, log2IV[s]) + "\n");
			}
		}
		double[] regCoeff = regressionLogVarianceVariedhn (numSkipRegression);
		//sb.append("  Slope of log2(var) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		//sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");
		sb.append("  C     = " + Math.exp(regCoeff[0]) + "\n");
		sb.append("  Slope of log2(IV) : beta  = " + -regCoeff[1] + "\n");
		sb.append("  delta = " + -regCoeff[2] + "\n");		
		//sb.append("  gamma = " + (-regCoeff[1])/(alpha - regCoeff[2]) + "\n");	
		//sb.append("  nu    = " + (-alpha * regCoeff[1])/(alpha - regCoeff[2]) + "\n\n");	
		sb.append("  Total CPU Time = " + cpuTime + "\n");
		sb.append("-----------------------------------------------------\n");		
		return sb.toString();
	}
	
	public String reportD(boolean details) {
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for density estimation: \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + numReplicates + "\n");
		sb.append(" Point sets: " + this.toString() + "\n\n");
		sb.append("RQMC variance \n");
		if (details) {
			sb.append("    n      log2(IV) \n");
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				sb.append(" " + size[s] + 
				          " " + PrintfFormat.f(7, 2, log2IV[s]) + "\n");
			}
		}
		double[] regCoeff = regressionLogVariance (numSkipRegression);
		
		sb.append("  Slope of log2(IV) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");
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
         	this.testVarianceRate (model, m, listDE.get(i), numEvalPoints, integVariance, ptSeries);
			sb.append (reportD (details));	
			makePlotsVar (numSets, m, (model.toString()).split(" ")[0], " ");
			//System.out.println("test");
		  }
	    }
		return sb.toString();
	}
	
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
         	//if (listDE.get(i)== new DEHistogram(listDE.get(i).getA(),listDE.get(i).getB()))
			sb.append (reportVarianceVariedhn (details));	
			makePlotsVar (numSets, m, (model.toString()).split(" ")[0], " ");
		  }
	    }
		return sb.toString();
	}
	
	public void makePlotsVar (int numSets, int m, String descModel, String descPoints) {
		// makeGraph();
		try {
			// String title = descModel + "; " + descPoints;
			// double [][] bidon = new double[0][0];
			// XYLineChart chart = new XYLineChart (title, "lg n", "lg MISE", bidon);
			// chart.init (title, "lg n", "lg MISE");
		    
			//for (int j = 3; j < numStats; j++)
				// chart.add(log2n, log2StatsMISE[j], statNames[j], " ");
			    chart.add(log2n, log2IV, " ", " ");
			FileWriter file = new FileWriter(descModel + "_" + descPoints + "_Var.tex");
			file.write(chart.toLatex(12, 8));
			file.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Performs an experiment (testVarianceRate) for each point set series in the given list,
	 * and returns a report as a string. 
	 * 
	 */
	
	public void fitPrintRegressionVariance (int order,ArrayList<DensityEstimator> listDE, boolean useBandwidth, 
			  int start, int number,
			  double range, double alpha, String method, StringBuffer sb) {
		double[][] regDataX = new double[listDE.size() * number][order];
		double[] regDataY = new double[listDE.size()  * number];
		double[] regDataYVar = new double[listDE.size()  * number];
		double[] coef;
		double logh;
		for (int j = 0; j < listDE.size() ; j++) { // For each m.
				logh = Num.log2(listDE.get(j).geth());  // For KDE.	
			for (int s = 0; s < number; s++) { // For each cardinality n
				regDataX[j * number + s][0] = log2n[start+s];
				regDataX[j  * number + s][1] = logh;
				if (order > 2) regDataX[j * number + s][2] = logh * log2n[start+s];
				regDataYVar[j * number + s] = log2VarS[j][start+s];
			}
		}
		/*coef = LeastSquares.calcCoefficients0(regDataX, regDataY);
		// System.out.println(" coef computed 1.\n");
		sb.append("  Regression coefficients for " + method + ".\n");
		sb.append("  C     = " + Math.exp(coef[0]) + "\n");
		sb.append("  beta  = " + -coef[1] + "\n");
		sb.append("  delta = " + -coef[2] + "\n");
		if (order > 2) sb.append("  inter = " + -coef[3] + "\n");
		sb.append("  gamma = " + (-coef[1])/(alpha - coef[2]) + "\n");	
		sb.append("  nu    = " + (-alpha * coef[1])/(alpha - coef[2]) + "\n\n");	*/

			coef = LeastSquares.calcCoefficients0(regDataX, regDataYVar);
			sb.append("  C for MISE     = " + Math.exp(coef[0]) + "\n");
			sb.append("  beta for MISE  = " + -coef[1] + "\n");
			sb.append("  delta for MISE = " + -coef[2] + "\n\n");
		}
	
	public String toString () {
		return theSets[0].toString();
	}
}