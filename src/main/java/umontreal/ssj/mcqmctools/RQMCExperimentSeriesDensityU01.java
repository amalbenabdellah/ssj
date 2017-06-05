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
import umontreal.ssj.stat.TallyStore;
import umontreal.ssj.stat.density.DEAveragedShiftedHistogram;
import umontreal.ssj.stat.density.DEHistogram;
import umontreal.ssj.stat.density.DEKernelDensity;
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

public class RQMCExperimentSeriesDensityU01 extends RQMCExperimentSeriesDensityKnown {


   /**
    * Constructor with a give series of RQMC point sets.
    *  @param theSets      the RQMC point sets
    */
	XYLineChart chart = new XYLineChart ();
	
double[] log2MISE = new double[numSets]; // log_2 of MISE
public double[][] meanD ;
public double[][] log2MISES;
public double[] log2KS ;
public double[] log2CVM ;
	
	//double[] log2ISB = new double[numSets]; // log_2 of h
   public RQMCExperimentSeriesDensityU01 (RQMCPointSet[] theSets) {
	   super(theSets);
	   log2MISE = new double[numSets]; // log_2 of the MISE
	   log2KS = new double[numSets];
	   log2CVM = new double[numSets];
   }
   

   /**
    * Performs an RQMC experiment with the given model for fixed h , with this series of RQMC point sets.  
    * For each set in the series, the variance, the integrated variance (IV), the kolmogorov smirnov test (KS), the Cramer von Mises (CvM) and   the mean integrated square error (MISE), its log in base 2.
    */
   
   
   public void testMISERate (MonteCarloModelDensityU01 model, int m,
			DensityEstimator DE, int numEvalPoints,  double[] MISE, double[] integVariance, double[] bias, RQMCPointSet [] theSets) {
	int n;
	
	Chrono timer = new Chrono();
	numReplicates = m;
	this.model = model;
	Tally statReps = new Tally();
	TallyStore statKS = new TallyStore();
	TallyStore statCVM = new TallyStore();
   if (displayExec) {
   	System.out.println("\n ============================================= ");
   	System.out.println("RQMC simulation for density estimation :  ");
   	System.out.println("Model: " + model.toString());
   	System.out.println(" Number of indep copies m  = " + m);
   	System.out.println(" Point sets: " + theSets[0].toString() + "\n");
    System.out.println("    n    log2Var log2KS log2CVM  log2(IV)  log2(MISE)   \n");    	
   }
  
	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		size[s] = n;
		double[][] data = new double[m][];		
		log2n[s] = Num.log2(n);
		RQMCExperiment.simulReplicatesRQMCSaveU01 (model, theSets[s], m, statReps, statKS, statCVM, data);
		integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model,  m, data, DE, numEvalPoints);
		MISE[s]=RQMCExperimentDensity.computeDensityMISE (model, m, data, DE, numEvalPoints);
		bias[s]=RQMCExperimentDensity.computeDensityBias (model,  m, data, DE, numEvalPoints);
		log2KS[s] = Num.log2(statKS.average());
		log2CVM[s] = Num.log2(statCVM.average());
		System.out.println("KS"+log2KS[s]+"CVM"+log2CVM[s]);
		log2Var[s] = Num.log2(statReps.variance());
	    log2MISE[s] = Num.log2(MISE[s]);
	    log2IV[s] = Num.log2(integVariance[s]);
	    if (displayExec) {

	    System.out.println( " " + size[s] + " " + PrintfFormat.f(10, 5, log2Var[s]) + " " + PrintfFormat.f(10, 5, log2KS[s]) + " " + PrintfFormat.f(10, 5, log2CVM[s]) +" " + PrintfFormat.f(10, 5, log2IV[s]) +
	          " " + PrintfFormat.f(7, 2, log2MISE[s]) + "\n");
	    }
	}	
	 
   cpuTime = timer.format();	   
}

   
   /**
    * Performs an RQMC experiment with the given model for varied h and n , with this series of RQMC point sets.  
   * For each set in the series, computes the bias, the variance, the integrated variance (IV), the integrated square bias (ISB), the kolmogorov smirnov test (KS), the Cramer von Mises (CvM) and   the mean integrated square error (MISE), its log in base 2.
    */

   public void testMISERateVariedhn (MonteCarloModelDensityU01 model, int m,
			DensityEstimator DE, int numEvalPoints,  double[] MISE, double[] integVariance, double[] bias, RQMCPointSet [] theSets) {
	int n;
	Tally statReps = new Tally();
	TallyStore statKS = new TallyStore();
	TallyStore statCVM = new TallyStore();
	Chrono timer = new Chrono();
	numReplicates = m;
	this.model = model;
  if (displayExec) {
  	System.out.println("\n ============================================= ");
  	System.out.println("RQMC simulation for density estimation :  ");
  	System.out.println("Model: " + model.toString());
  	System.out.println(" Number of indep copies m  = " + m);
  	System.out.println(" Point sets: " + theSets[0].toString() + "\n");
  	 System.out.println("    n    log2Var log2KS log2CVM  log2(IV) log2(ISB)  log2(MISE)   \n");   	    	
  }

  double  l=1;

	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		size[s] = n;
		double[][] data = new double[m][];	
		log2n[s] = Num.log2(n);
		RQMCExperiment.simulReplicatesRQMCSaveU01 (model, theSets[s], m, statReps, statKS, statCVM, data);			
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
			l=l+0.2;
		}
		
		
		
		integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model,  m, data, DE, numEvalPoints);
		MISE[s]=RQMCExperimentDensity.computeDensityMISE (model, m, data, DE, numEvalPoints);
		bias[s]=RQMCExperimentDensity.computeDensityBias (model,  m, data, DE, numEvalPoints);
		log2Var[s] = Num.log2(statReps.variance());
		log2KS[s] = Num.log2(statKS.average());
		log2CVM[s] = Num.log2(statCVM.average());
	    log2MISE[s] = Num.log2(MISE[s]);
	    log2IV[s] = Num.log2(integVariance[s]);
	    log2ISB[s]= Num.log2(bias[s]);
	    if (displayExec) {
	    	System.out.println( " " + size[s] + " " + PrintfFormat.f(10, 5, log2Var[s]) + " " + PrintfFormat.f(10, 5, log2KS[s]) + " " + PrintfFormat.f(10, 5, log2CVM[s]) +" " + PrintfFormat.f(10, 5, log2IV[s]) + " " + PrintfFormat.f(10, 5, log2ISB[s]) +
	  	          " " + PrintfFormat.f(7, 2, log2MISE[s]) + "\n");
	    }
	}	
	 
  cpuTime = timer.format();	   
}
   
   /**
    * Performs an RQMC experiment with the given model for h chosen optimally in function of n , with this series of RQMC point sets.  
   * For each set in the series, computes the bias, the variance, the integrated variance (IV), the integrated square bias (ISB), the kolmogorov smirnov test (KS), the Cramer von Mises (CvM) and   the mean integrated square error (MISE), its log in base 2.
   */
   public void testMISERateOptimal (MonteCarloModelDensityU01 model, int m,
			DensityEstimator DE, int numEvalPoints,  double[] MISE, double[] integVariance, double[] bias, RQMCPointSet [] theSets) {
	int n;
	int r = 4;
	Tally statReps = new Tally();
	TallyStore statKS = new TallyStore();
	TallyStore statCVM = new TallyStore();
	Chrono timer = new Chrono();
	numReplicates = m;
	this.model = model;
  if (displayExec) {
  	System.out.println("\n ============================================= ");
  	System.out.println("RQMC simulation for density estimation:  ");
  	System.out.println("Model: " + model.toString());
  	System.out.println(" Number of indep copies m  = " + m);
  	System.out.println(" Point sets: " + theSets[0].toString() + "\n");
  	 System.out.println("    n    log2Var log2KS log2CVM  log2(IV) log2(ISB)  log2(MISE)   \n");       	
  }
  
  testMISERateVariedhn(model, m, DE, numEvalPoints,MISE,integVariance,bias, theSets); 
  
  double[] regCoefbias= slope (log2h, log2ISB, numSkipRegression);
	double alpha=regCoefbias[1];
	System.out.println("alpha"+regCoefbias[1] );

	double[] regCoeff = slope( log2n, log2h, log2IV, numSkipRegression);
	double C=Math.exp(regCoeff[0]);
	double B=Math.exp(regCoefbias[0]);
	
	
	double  beta  =  -regCoeff[1] ;
	double delta =   -regCoeff[2];	
	double gamma = beta/(alpha+delta);	
	double kappa= Math.pow((C*delta/B*alpha),1/(alpha+delta));
	System.out.println("gamma"+gamma );
 
  
 
	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		size[s] = n;
		double[][] data = new double[m][];		
		log2n[s] = Num.log2(n);
		DE.seth(kappa*Math.pow(n, -gamma));
		//DE.seth(1/Math.pow(8, r)*Math.pow(n, -gamma));
		RQMCExperiment.simulReplicatesRQMCSaveU01 (model, theSets[s], m, statReps, statKS, statCVM, data);	
		integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model,  m, data, DE, numEvalPoints);
		MISE[s]=RQMCExperimentDensity.computeDensityMISE (model, m, data, DE, numEvalPoints);
		bias[s]=RQMCExperimentDensity.computeDensityBias (model,  m, data, DE, numEvalPoints);
		log2Var[s] = Num.log2(statReps.variance());
		log2KS[s] = Num.log2(statKS.average());
		log2CVM[s] = Num.log2(statCVM.average());
		System.out.println("MI"+MISE[s] );
		//mean[s] = statReps.average();
	    log2MISE[s] = Num.log2(MISE[s]);
	    log2IV[s] = Num.log2(integVariance[s]);
	    log2ISB[s]= Num.log2(bias[s]);
	    if (displayExec) {
	    	System.out.println( " " + size[s] + " " + PrintfFormat.f(10, 5, log2Var[s]) + " " + PrintfFormat.f(10, 5, log2KS[s]) + " " + PrintfFormat.f(10, 5, log2CVM[s]) +" " + PrintfFormat.f(10, 5, log2IV[s]) + " " + PrintfFormat.f(10, 5, log2ISB[s]) +
		  	          " " + PrintfFormat.f(7, 2, log2MISE[s]) + "\n");
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
    * Produces and returns a report on the last experiment for varied h and n.
    * @param details  If true, gives values ( logVar, logIV, logMISE, logKS, logCvM...) for each n.
    * @param densityestimator, the name of the statistic used 
    * @return  Report as a string.
    */
	
	
	
	public String reportMISEVAriedhn (boolean details , String densityestimator) {
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for density estimation: \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + numReplicates + "\n");
		sb.append(" Density Esimator  = " + densityestimator + "\n");
		sb.append(" Point sets: " + this.toString() + "\n\n");		
		sb.append("RQMC Mean Integrated Square Error (MISE) \n");
		if (details) {
			sb.append("    n    log2Var log2KS log2CVM  log2(IV) log2(ISB)  log2(MISE)   \n");  
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				sb.append( " " + size[s] + " " + PrintfFormat.f(10, 5, log2Var[s]) + " " + PrintfFormat.f(10, 5, log2KS[s]) + " " + PrintfFormat.f(10, 5, log2CVM[s]) +" " + PrintfFormat.f(10, 5, log2IV[s]) + " " + PrintfFormat.f(10, 5, log2ISB[s]) +
			  	          " " + PrintfFormat.f(7, 2, log2MISE[s]) + "\n");
			}
		}
		double[] regCoefbias= slope (log2h, log2ISB, numSkipRegression);
		double alpha=regCoefbias[1];
		double[] regCoeff = slope (log2n, log2h, log2MISE, numSkipRegression);
		sb.append("  beta for MISE = " + -regCoeff[1] + "\n");
		double delta = alpha - regCoeff[2];
		sb.append("  delta for MISE = " + delta + "\n");			
		sb.append("  gamma = " + (-regCoeff[1])/(alpha - regCoeff[2]) + "\n");	
		sb.append("  nu    = " + (-alpha * regCoeff[1])/(alpha - regCoeff[2]) + "\n\n");
		
		regCoeff = slope(log2n, log2IV, numSkipRegression);
		sb.append("  Slope of log2(IV) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");
		regCoeff = slope(log2n, log2ISB, numSkipRegression);
		sb.append("  Slope of log2(ISB) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		regCoeff = slope(log2n, log2Var, numSkipRegression);
		sb.append("  Slope of log2(VAr) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		//sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");
		regCoeff = slope(log2n, log2KS, numSkipRegression);
		sb.append("  Slope of log2(KS) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		//sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");
		regCoeff = slope(log2n, log2CVM, numSkipRegression);
		sb.append("  Slope of log2(CVM) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		//sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");
		sb.append("  Total CPU Time = " + cpuTime + "\n");
		sb.append("-----------------------------------------------------\n");		
		return sb.toString();
	}
	
	/**
	    * Produces and returns a report on the last experiment for fixed h.
	    * @param details  If true, gives values (logVar, logIV, logMISE, logKS, logCvM...) for each n.
	    * @param densityestimator, the name of the statistic used 
	    * @return  Report as a string.
	    */
		
	public String reportMISE (boolean details, String densityestimator) {
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for density estimation : \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + numReplicates + "\n");
		sb.append(" Density Esimator  = " + densityestimator + "\n");
		sb.append(" Point sets: " + this.toString() + "\n\n");
		sb.append("RQMC Mean Integrated Square Error (MISE) \n");
		if (details) {
			sb.append("    n    log2Var log2KS log2CVM  log2(IV)   log2(MISE)   \n");  
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				sb.append( " " + size[s] + " " + PrintfFormat.f(10, 5, log2Var[s]) + " " + PrintfFormat.f(10, 5, log2KS[s]) + " " + PrintfFormat.f(10, 5, log2CVM[s]) +" " + PrintfFormat.f(10, 5, log2IV[s])  +
			  	          " " + PrintfFormat.f(7, 2, log2MISE[s]) + "\n");
			}
		}
		double[] regCoeff = slope(log2n, log2MISE, numSkipRegression);
		sb.append("  Slope of log2(MISE) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");
		regCoeff = slope(log2n, log2IV, numSkipRegression);
		sb.append("  Slope of log2(IV) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");
		regCoeff = slope(log2n, log2Var, numSkipRegression);
		sb.append("  Slope of log2(Var) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		regCoeff = slope(log2n, log2KS, numSkipRegression);
		sb.append("  Slope of log2(KS) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		regCoeff = slope(log2n, log2CVM, numSkipRegression);
		sb.append("  Slope of log2(CVM) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");		
		sb.append("  Total CPU Time = " + cpuTime + "\n");
		sb.append("-----------------------------------------------------\n");		
		return sb.toString();
	}
	
	
	/**
	 * Performs an experiment (testVarianceRate)  with the given model for each point set series in the given list,
	 * and returns a report as a string. 
	 * 
	 */
	
	
	public String TestRQMCManyPointTypes (MonteCarloModelDensityKnown model, 
			ArrayList<RQMCPointSet[]> list, int m,
			ArrayList<DensityEstimator> listDE, int numEvalPoints, 
            boolean details) {
		StringBuffer sb = new StringBuffer("");
		numReplicates = m;	
		double[] MISE= new double[numSets];   // Will contain the IV estimates, for each n.
		double[] IntegVariance= new double[numSets];
		double[] bias= new double[numSets];
		for(int i=0; i < listDE.size(); i++) {
		  for (RQMCPointSet[] ptSeries : list) {	
         	testMISERate (model, m, listDE.get(i), numEvalPoints, MISE, IntegVariance,bias, ptSeries);         	
  			sb.append (reportMISE (details, listDE.get(i).toString()));	
  			makePlotsMISE (numSets, m, (model.toString()).split(" ")[0], listDE.get(i).toString(), ptSeries.toString());
			
		  }
	    }
		return sb.toString();
	}
	public String TestRQMCManyPointTypesOptimal (MonteCarloModelDensityKnown model, 
			ArrayList<RQMCPointSet[]> list, int m,
			ArrayList<DensityEstimator> listDE, int numEvalPoints, 
            boolean details) {
		StringBuffer sb = new StringBuffer("");
		numReplicates = m;	
		double[] MISE= new double[numSets];   // Will contain the IV estimates, for each n.
		double[] IntegVariance= new double[numSets];
		double[] bias= new double[numSets];
		for(int i=0; i < listDE.size(); i++) {
		  for (RQMCPointSet[] ptSeries : list) {			
         	testMISERateOptimal (model, m, listDE.get(i), numEvalPoints, MISE, IntegVariance,bias, ptSeries);
  			sb.append (reportMISE (details, listDE.get(i).toString()));	
  			makePlotsMISE (numSets, m, (model.toString()).split(" ")[0], listDE.get(i).toString(), ptSeries.toString());
			
		  }
	    }
		return sb.toString();
	}
	
	public String TestRQMCManyPointTypesVariedhn (MonteCarloModelDensityKnown model, 
			ArrayList<RQMCPointSet[]> list, int m,
			ArrayList<DensityEstimator> listDE, int numEvalPoints, 
            boolean details) {
		StringBuffer sb = new StringBuffer("");
		numReplicates = m;	
		double[] MISE= new double[numSets]; //Will contain the IV estimates, for each n.
		double[] integVariance= new double[numSets];
		double[] bias= new double[numSets];
		for(int i=0; i < listDE.size(); i++) {
		  for (RQMCPointSet[] ptSeries : list) {
         	testMISERateVariedhn (model, m, listDE.get(i), numEvalPoints, MISE, integVariance,bias, ptSeries);         	
			sb.append ( reportMISEVAriedhn (details, listDE.get(i).toString()));	
			makePlotsMISE (numSets, m, (model.toString()).split(" ")[0], listDE.get(i).toString(), ptSeries.toString());
		  }
	    }
		return sb.toString();
	}
	
	
	/** Method that produces Latex plots for the experiment for MISE
	 */
	
	public void makePlotsMISE (int numSets, int m, String descModel, String densityEstimator, String descPoints) {

		try {
			
			    chartMISE.add(log2n, log2MISE, " ", " ");
			FileWriter file = new FileWriter(descModel + "_" +   densityEstimator+ "_" + descPoints + "_MISE.tex");
			file.write(chartMISE.toLatex(12, 8));
			file.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	

	public String toString () {
		return theSets[0].toString();
	}
}