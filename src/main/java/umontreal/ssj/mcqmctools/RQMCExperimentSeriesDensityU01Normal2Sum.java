/*package umontreal.ssj.mcqmctools;

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
 


import umontreal.ssj.charts.XYLineChart;
import umontreal.ssj.functionfit.LeastSquares;
import umontreal.ssj.hups.*;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stat.TallyStore;
import umontreal.ssj.stat.density.DEAveragedShiftedHistogram;
import umontreal.ssj.stat.density.DEAveragedShiftedHistogramWeight;
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


*//**
 * This class offers facilities to perform experiments to study the convergence
 * of the variance when estimating a mean (expectation) with a series of RQMC 
 * point sets usually of the same type, but different size @f$n@f.
 * The series of RQMC point sets of different sizes can be passed in an array 
 * to the constructor.   The method @f$testVarianceRate@f$ performs an experiment 
 * with a given model and the series of point sets.  One can recover the average,
 * variance, their logs in base 2, etc., in arrays, as well as the estimated 
 * linear regression of log(variance) as a function of log(n). 
 *//*

public class RQMCExperimentSeriesDensityU01 extends RQMCExperimentSeriesDensityKnown {


   *//**
    * Constructor with a give series of RQMC point sets.
    *  @param theSets      the RQMC point sets
    *//*
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
   

   *//**
    * Performs an RQMC experiment with the given model for fixed h , with this series of RQMC point sets.  
    * For each set in the series, the variance, the integrated variance (IV), the kolmogorov smirnov test (KS), the Cramer von Mises (CvM) and   the mean integrated square error (MISE), its log in base 2.
    *//*
   
   
   public void teststatistics (MonteCarloModelDensityU01 model, int m, RQMCPointSet [] theSets) {
	   int n;
	   Tally statReps = new Tally();
	   Chrono timer = new Chrono();
	   System.out.println("    n    log2(KS) log2(CVM)  log2(Var)   \n");  
	   for (int s = 0; s < numSets; s++) { // For each cardinality n
			n = theSets[s].getNumPoints();
			TallyStore statKS = new TallyStore(n);
			TallyStore statCVM = new TallyStore(n);
			size[s] = n;
			double[][] data = new double[m][];		
			log2n[s] = Num.log2(n);
			RQMCExperiment.simulReplicatesRQMCSaveU01 (model, theSets[s], m, statReps, statKS, statCVM, data);
			log2KS[s] = Num.log2(statKS.average());
			log2CVM[s] = Num.log2(statCVM.average());
			log2Var[s] = Num.log2(statReps.variance());
			
		    if (displayExec) {
		    	

		    	System.out.println( " " + size[s] + " "  + PrintfFormat.f(10, 5, log2KS[s]) + " " + PrintfFormat.f(10, 5, log2CVM[s]) +" " + PrintfFormat.f(10, 5, log2Var[s]) +"\n");
		    }
		}	
		 
	   cpuTime = timer.format();	
	   
	   
	   
	   
	   
   }
   
   public void testMISERateU01 (MonteCarloModelDensityU01 model, int m,
			DensityEstimator DE, int numEvalPoints,  double[] MISE, double[] integVariance, double[] bias, RQMCPointSet [] theSets) {
	int n;
	
	Chrono timer = new Chrono();
	numReplicates = m;
	this.model = model;
	Tally statReps = new Tally();
	  double  t=0.05;
	
   if (displayExec) {
   	System.out.println("\n ============================================= ");
   	System.out.println("RQMC simulation for density estimation with fixed h:  ");
   	System.out.println("Model: " + model.toString());
   	System.out.println(" Number of indep copies m  = " + m);
   	System.out.println(" Point sets: " + theSets[0].toString() + "\n");
   // System.out.println("    n    log2Var log2KS log2CVM  log2(IV)  log2(MISE)   \n");   
    System.out.println("    n    log2(IV)  log2(MISE)   \n");   
   }
  
	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		TallyStore statKS = new TallyStore(n);
		TallyStore statCVM = new TallyStore(n);
		size[s] = n;
		double[][] data = new double[m][];		
		log2n[s] = Num.log2(n);
		log2h[s] = Num.log2(t);
		t+=0.05;
		//RQMCExperiment.simulReplicatesRQMCSaveU01 (model, theSets[s], m, statReps, statKS, statCVM, data);
		RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statReps, data);
		integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model,  m, data, DE, numEvalPoints);
		MISE[s]=RQMCExperimentDensity.computeDensityMISE (model, m, data, DE, numEvalPoints);
		bias[s]=RQMCExperimentDensity.computeDensityBias (model,  m, data, DE, numEvalPoints);
		//log2KS[s] = Num.log2(statKS.average());
		//log2CVM[s] = Num.log2(statCVM.average());
		//log2Var[s] = Num.log2(statReps.variance());
	    log2MISE[s] = Num.log2(MISE[s]);
	    log2IV[s] = Num.log2(integVariance[s]);
	    log2ISB[s] = Num.log2(bias[s]);
	    if (displayExec) {

	    	System.out.println( " " + size[s] + " " + PrintfFormat.f(10, 5, log2Var[s]) + " " + PrintfFormat.f(10, 5, log2KS[s]) + " " + PrintfFormat.f(10, 5, log2CVM[s]) +" " + PrintfFormat.f(10, 5, log2IV[s]) +
		  	          " " + PrintfFormat.f(10, 5, log2MISE[s]) + "\n");
	    	System.out.println( " " + size[s]  +" " + PrintfFormat.f(10, 5, log2IV[s])  +
		  	          " " + PrintfFormat.f(10, 5, log2MISE[s]) + "\n");
	    }
	}	
	 
   cpuTime = timer.format();	   
}

   
   *//**
    * Performs an RQMC experiment with the given model for varied h and n , with this series of RQMC point sets.  
   * For each set in the series, computes the bias, the variance, the integrated variance (IV), the integrated square bias (ISB), the kolmogorov smirnov test (KS), the Cramer von Mises (CvM) and   the mean integrated square error (MISE), its log in base 2.
    *//*

   public void testMISERateVariedhnU01 (MonteCarloModelDensityU01 model, int m,
			DensityEstimator DE, int numEvalPoints,  double[] MISE, double[] integVariance, double[] bias, RQMCPointSet [] theSets) {
	int n;
	Tally statReps = new Tally();
	
	Chrono timer = new Chrono();
	numReplicates = m;
	this.model = model;
  if (displayExec) {
  	System.out.println("\n ============================================= ");
  	System.out.println("RQMC simulation for density estimation with varied h and n :  ");
  	System.out.println("Model: " + model.toString());
  	System.out.println(" Number of indep copies m  = " + m);
  	System.out.println(" Point sets: " + theSets[0].toString() + "\n");
  	// System.out.println("    n    log2Var log2KS log2CVM  log2(IV) log2(ISB)  log2(MISE)   \n");   	
  	System.out.println("    n    log2(IV) log2(ISB)  log2(MISE)   \n");   	   
  }

  double  l=0.01, t=0.05;

	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		TallyStore statKS = new TallyStore(n);
		TallyStore statCVM = new TallyStore(n);
		size[s] = n;
		double[][] data = new double[m][];	
		log2n[s] = Num.log2(n);
		//RQMCExperiment.simulReplicatesRQMCSaveU01 (model, theSets[s], m, statReps, statKS, statCVM, data);
				RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statReps, data);	
		if ( DE == new DEHistogram(DE.getA(),DE.getB()) ){
		    log2h[s] = Num.log2((DE.getB()-DE.getA())/Math.pow(4, l));
		    DE.seth((DE.getB()-DE.getA())/Math.pow(4, l));
		    l++;
			}
		
		else if ( DE == new DEAveragedShiftedHistogram(DE.getA(),DE.getB()) || DE == new DEAveragedShiftedHistogramWeight(DE.getA(),DE.getB())){
			log2h[s] = Num.log2((DE.getB()-DE.getA())/(32*Math.pow(4, l)));
		    DE.seth((DE.getB()-DE.getA())/(32*Math.pow(4, l)));	
		    l++;
			
		}
		
		if ( DE == new DEHistogram(DE.getA(),DE.getB()) ||DE == new DEAveragedShiftedHistogram(DE.getA(),DE.getB()) || DE == new DEAveragedShiftedHistogramWeight(DE.getA(),DE.getB()) ){
		    log2h[s] = Num.log2((DE.getB()-DE.getA())/Math.pow(4, l));
		    DE.seth((DE.getB()-DE.getA())/Math.pow(4, l));
		    l++;
			}
		else {			
			log2h[s] =- Num.log2(Math.pow(4, l));
			 DE.seth(1/Math.pow(4, l));
			l=l+0.2;
		}
		if (DE == new DEKernelDensity(DE.getA(),DE.getB())){
		
		log2h[s] = Num.log2(t);
		DE.seth(t);
		t+=0.05;	}
		
		else
			
		{
			log2h[s] = Num.log2((model.getMax()-model.getMin())/t);
			DE.seth((model.getMax()-model.getMin())/t);
			t+=0.05;	
			
		}
		if (DE == new DEKernelDensity(DE.getA(),DE.getB())){
		log2h[s] = Num.log2(t);
		DE.seth(t);
		t+=0.05;}
		else{
		log2h[s] = Num.log2(l);
		DE.seth(l);
		l+=0.02;}
		
		log2h[s] = Num.log2(t);
		DE.seth(t);
		t+=0.05;
				
				log2h[s] = Num.log2(t);
				DE.seth(t);
				t+=0.05;
		
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
	  	          " " + PrintfFormat.f(10, 5, log2MISE[s]) + "\n");
	    	System.out.println( " " + size[s]  +" " + PrintfFormat.f(10, 5, log2IV[s]) + " " + PrintfFormat.f(10, 5, log2ISB[s]) +
		  	          " " + PrintfFormat.f(10, 5, log2MISE[s]) + "\n");
	    }
	}	
	 
  cpuTime = timer.format();	   
}
   
   *//**
    * Performs an RQMC experiment with the given model for h chosen optimally in function of n , with this series of RQMC point sets.  
   * For each set in the series, computes the bias, the variance, the integrated variance (IV), the integrated square bias (ISB), the kolmogorov smirnov test (KS), the Cramer von Mises (CvM) and   the mean integrated square error (MISE), its log in base 2.
   *//*
   public void testMISERateOptimalU01 (MonteCarloModelDensityU01 model, int m,
			DensityEstimator DE, int numEvalPoints,  double[] MISE, double[] integVariance, double[] bias, RQMCPointSet [] theSets) {
	int n;
	int r = 4;
	Tally statReps = new Tally();
	
	Chrono timer = new Chrono();
	numReplicates = m;
	this.model = model;
  if (displayExec) {
  	System.out.println("\n ============================================= ");
  	System.out.println("RQMC simulation for density estimation with optimal h:  ");
  	System.out.println("Model: " + model.toString());
  	System.out.println(" Number of indep copies m  = " + m);
  	System.out.println(" Point sets: " + theSets[0].toString() + "\n");
  	 System.out.println("    n    log2(IV) log2(ISB)  log2(MISE)   \n");  
  	 
  }
  
  testMISERateU01(model, m, DE, numEvalPoints,MISE,integVariance,bias, theSets); 
  
  
  
  double[] regCoefbias= slope (log2h, log2ISB, numSkipRegression);
	double alpha=regCoefbias[1];
	System.out.println("alpha"+regCoefbias[1] );
  double alpha;
	 if(DE==new DEHistogram(DE.getA(),DE.getB()))
		 alpha = 2;
	 else
		 alpha = 4;

	double[] regCoeff = slope( log2n, log2h, log2IV, numSkipRegression);
	double C=Math.exp(regCoeff[0]);
	double B=Math.exp(regCoefbias[0]);
	
	
	double  beta  =  -regCoeff[1] ;
	double delta =   -regCoeff[2];	
	double gamma = beta/(alpha+delta);	
	double kappa= Math.pow((C*delta/B*alpha),1/(alpha+delta));
	System.out.println("gamma"+gamma );
	System.out.println("kappa"+kappa );
 
  
 
	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		size[s] = n;
		TallyStore statKS = new TallyStore(n);
		TallyStore statCVM = new TallyStore(n);
		double[][] data = new double[m][];		
		log2n[s] = Num.log2(n);
		DE.seth(kappa*Math.pow(n, -gamma));
		//DE.seth(1/Math.pow(8, r)*Math.pow(n, -gamma));
		//RQMCExperiment.simulReplicatesRQMCSaveU01 (model, theSets[s], m, statReps, statKS, statCVM, data);
				RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statReps, data);
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
	    	System.out.println( " " + size[s]  +" " + PrintfFormat.f(10, 5, log2IV[s]) + " " + PrintfFormat.f(10, 5, log2ISB[s]) +
		  	          " " + PrintfFormat.f(10, 5, log2MISE[s]) + "\n");
	    }
	}	
	 
  cpuTime = timer.format();	   
}

   *//**
    * Returns the vector of log_2(n).
    *//*
   public double[] getLog2n() {
      return log2n;
   }
   
   *//**
    * Produces and returns a report on the last experiment for varied h and n.
    * @param details  If true, gives values ( logVar, logIV, logMISE, logKS, logCvM...) for each n.
    * @param densityestimator, the name of the statistic used 
    * @return  Report as a string.
    *//*
	
	
	
	public String reportMISEVAriedhnU01 (boolean details , String densityestimator) {
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for density estimation: \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + numReplicates + "\n");
		sb.append(" Density Esimator  = " + densityestimator + "\n");
		sb.append(" Point sets: " + this.toString() + "\n\n");		
		sb.append("RQMC Mean Integrated Square Error (MISE) \n");
		if (details) {
			sb.append("    n      log2(IV) log2(ISB)  log2(MISE)   \n");  
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				sb.append( " " + size[s] + " " +  PrintfFormat.f(10, 5, log2IV[s]) + " " + PrintfFormat.f(10, 5, log2ISB[s]) +
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
	
	*//**
	    * Produces and returns a report on the last experiment for fixed h.
	    * @param details  If true, gives values (logVar, logIV, logMISE, logKS, logCvM...) for each n.
	    * @param densityestimator, the name of the statistic used 
	    * @return  Report as a string.
	    *//*
		
	public String reportMISEU01 (boolean details, String densityestimator) {
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for density estimation : \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + numReplicates + "\n");
		sb.append(" Density Esimator  = " + densityestimator + "\n");
		sb.append(" Point sets: " + this.toString() + "\n\n");
		sb.append("RQMC Mean Integrated Square Error (MISE) \n");
		if (details) {
			sb.append("    n     log2(IV)  log2(ISB) log2(MISE)   \n");  
			for (int s = 0; s < numSets; s++) { // For each cardinality n
					sb.append( " " + size[s] + " " + PrintfFormat.f(10, 5, log2IV[s]) + " " + PrintfFormat.f(10, 5, log2ISB[s]) +
				  	          " " + PrintfFormat.f(7, 2, log2MISE[s]) + "\n");
				}
			}
		
		double[] regCoeff = slope(log2n, log2MISE, numSkipRegression);
		sb.append("  Slope of log2(MISE) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");
		regCoeff = slope(log2n, log2IV, numSkipRegression);
		sb.append("  Slope of log2(IV) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");
		regCoeff = slope(log2n, log2ISB, numSkipRegression);
		sb.append("  Slope of log2(ISB) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
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
	
	
	*//**
	 * Performs an experiment (testVarianceRate)  with the given model for each point set series in the given list,
	 * and returns a report as a string. 
	 * 
	 *//*
	
	
	public String TestRQMCManyPointTypes (MonteCarloModelDensityU01 model, 
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
         	testMISERateU01 (model, m, listDE.get(i), numEvalPoints, MISE, IntegVariance,bias, ptSeries);         	
  			sb.append (reportMISEU01 (details, listDE.get(i).toString()));	
  			makePlotsMISE (numSets, m, (model.toString()).split(" ")[0], listDE.get(i).toString(), ptSeries.toString());
			
		  }
	    }
		for (RQMCPointSet[] ptSeries : list) {
			teststatistics (model,m, ptSeries);
			double[] regCoeff = slope(log2n, log2Var, numSkipRegression);
			sb.append("  Slope of log2(Var) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
			regCoeff = slope(log2n, log2KS, numSkipRegression);
			sb.append("  Slope of log2(KS) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
			regCoeff = slope(log2n, log2CVM, numSkipRegression);
			sb.append("  Slope of log2(CVM) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
			}
		return sb.toString();
	}
	public String TestRQMCManyPointTypesOptimal (MonteCarloModelDensityU01 model, 
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
         	testMISERateOptimalU01 (model, m, listDE.get(i), numEvalPoints, MISE, IntegVariance,bias, ptSeries);
  			sb.append (reportMISEU01 (details, listDE.get(i).toString()));	
  			makePlotsMISE (numSets, m, (model.toString()).split(" ")[0], listDE.get(i).toString(), ptSeries.toString());
			
		  }
	    }
		for (RQMCPointSet[] ptSeries : list) {
			teststatistics (model,m, ptSeries);
			double[] regCoeff = slope(log2n, log2Var, numSkipRegression);
			sb.append("  Slope of log2(Var) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
			regCoeff = slope(log2n, log2KS, numSkipRegression);
			sb.append("  Slope of log2(KS) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
			regCoeff = slope(log2n, log2CVM, numSkipRegression);
			sb.append("  Slope of log2(CVM) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
			}
		return sb.toString();
	}
	
	public String TestRQMCManyPointTypesVariedhn (MonteCarloModelDensityU01 model, 
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
         	testMISERateVariedhnU01 (model, m, listDE.get(i), numEvalPoints, MISE, integVariance,bias, ptSeries);         	
			sb.append ( reportMISEVAriedhnU01 (details, listDE.get(i).toString()));				
			makePlotsMISE (numSets, m, (model.toString()).split(" ")[0], listDE.get(i).toString(), ptSeries.toString());
		  }
	    }
		for (RQMCPointSet[] ptSeries : list) {
		teststatistics (model,m, ptSeries);
		double[] regCoeff = slope(log2n, log2Var, numSkipRegression);
		sb.append("  Slope of log2(Var) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		regCoeff = slope(log2n, log2KS, numSkipRegression);
		sb.append("  Slope of log2(KS) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		regCoeff = slope(log2n, log2CVM, numSkipRegression);
		sb.append("  Slope of log2(CVM) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		}
		return sb.toString();
	}
	
	
	*//** Method that produces Latex plots for the experiment for MISE
	 *//*
	
	public void makePlotsMISE (int numSets, int m, String descModel, String densityEstimator, String descPoints) {

		try {
			
			    chartMISE.add(log2n, log2MISE, " ", " ");
			FileWriter file = new FileWriter(descModel + "_" +   densityEstimator+ "_" + descPoints + "_MISEU01.tex");
			file.write(chartMISE.toLatex(12, 8));
			file.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	

	public String toString () {
		return theSets[0].toString();
	}
}*/


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
import umontreal.ssj.stat.density.DEAveragedShiftedHistogramWeight;
import umontreal.ssj.stat.density.DEGaussianKDEBotev;
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

public class RQMCExperimentSeriesDensityU01Normal2Sum extends RQMCExperimentSeriesDensityKnown {


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
	
   public RQMCExperimentSeriesDensityU01Normal2Sum (RQMCPointSet[] theSets) {
	   super(theSets);
	   log2MISE = new double[numSets]; // log_2 of the MISE
	   log2KS = new double[numSets];
	   log2CVM = new double[numSets];
   }
   

   /**
    * Performs an RQMC experiment with the given model for fixed h , with this series of RQMC point sets.  
    * For each set in the series, the variance, the integrated variance (IV), the kolmogorov smirnov test (KS), the Cramer von Mises (CvM) and   the mean integrated square error (MISE), its log in base 2.
    */
   
   public void teststatistics (MonteCarloModelDensityU01 model, int m, RQMCPointSet [] theSets) {
	   int n;
	   Tally statReps = new Tally();
	   Chrono timer = new Chrono();
	   System.out.println("    n    log2(KS) log2(CVM)  log2(Var)   \n");  
	   for (int s = 0; s < numSets; s++) { // For each cardinality n
			n = theSets[s].getNumPoints();
			TallyStore statKS = new TallyStore(n);
			TallyStore statCVM = new TallyStore(n);
			size[s] = n;
			double[][] data = new double[m][];		
			log2n[s] = Num.log2(n);
			RQMCExperiment.simulReplicatesRQMCSaveU01 (model, theSets[s], m, statReps, statKS, statCVM, data);
			log2KS[s] = Num.log2(statKS.average());
			log2CVM[s] = Num.log2(statCVM.average());
			log2Var[s] = Num.log2(statReps.variance());
			
		    if (displayExec) {
		    	

		    	System.out.println( " " + size[s] + " "  + PrintfFormat.f(10, 5, log2KS[s]) + " " + PrintfFormat.f(10, 5, log2CVM[s]) +" " + PrintfFormat.f(10, 5, log2Var[s]) +"\n");
		    }
		}	
		 
	   cpuTime = timer.format();	
	   
	   
	   
	   
	   
   }
   
   public void testMISERateU01 (MonteCarloModelDensityU01 model, int m,
			DensityEstimator DE, int numEvalPoints,  double[] MISE, double[] integVariance, double[] bias, RQMCPointSet [] theSets) {
	int n;
	
	Chrono timer = new Chrono();
	numReplicates = m;
	this.model = model;
	Tally statReps = new Tally();
	
   if (displayExec) {
   	System.out.println("\n ============================================= ");
   	System.out.println("RQMC simulation for density estimation with fixed h:  ");
   	System.out.println("Model: " + model.toString());
   	System.out.println(" Number of indep copies m  = " + m);
   	System.out.println(" Point sets: " + theSets[0].toString() + "\n");
   	System.out.println("    n    log2(IV)   log2(MISE)   \n"); 
   }
  
	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		size[s] = n;
		double[][] data = new double[m][];		
		log2n[s] = Num.log2(n);
		RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statReps, data);
		integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model,  m, data, DE, numEvalPoints);
		MISE[s]=RQMCExperimentDensity.computeDensityMISE (model, m, data, DE, numEvalPoints);
		bias[s]=RQMCExperimentDensity.computeDensityBias (model,  m, data, DE, numEvalPoints);
	    log2MISE[s] = Num.log2(MISE[s]);
	    log2IV[s] = Num.log2(integVariance[s]);
	    if (displayExec) {
	    	System.out.println( " " + size[s] + " " + PrintfFormat.f(10, 5, log2IV[s]) +
		  	          " " + PrintfFormat.f(10, 5, log2MISE[s]) + "\n");
	    }
	}	
	 
   cpuTime = timer.format();	   
}

   
   /**
    * Performs an RQMC experiment with the given model for varied h and n , with this series of RQMC point sets.  
   * For each set in the series, computes the bias, the variance, the integrated variance (IV), the integrated square bias (ISB), the kolmogorov smirnov test (KS), the Cramer von Mises (CvM) and   the mean integrated square error (MISE), its log in base 2.
    */

   public void testMISERateVariedhnU01 (MonteCarloModelDensityU01 model, int m,
			DensityEstimator DE, int numEvalPoints,  double[] MISE, double[] integVariance, double[] bias, RQMCPointSet [] theSets, String point) {
	int n;
	Tally statReps = new Tally();
	
	Chrono timer = new Chrono();
	numReplicates = m;
	this.model = model;
  if (displayExec) {
  	System.out.println("\n ============================================= ");
  	System.out.println("RQMC simulation for density estimation with varied h and n :  ");
  	System.out.println("Model: " + model.toString());
  	System.out.println(" Number of indep copies m  = " + m);
  	System.out.println(" Point sets: " + point + "\n");
  	System.out.println("    n    log2(IV) log2(ISB)  log2(MISE)   \n"); 
  }

 // double t=0.009, l=0.055;
  
double r=6, d=2;
double l=0.0015;


	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		size[s] = n;
		double[][] data = new double[m][];	
		log2n[s] = Num.log2(n);		
		
		/*if ( DE == new DEHistogram(DE.getA(),DE.getB()) ||DE == new DEAveragedShiftedHistogram(DE.getA(),DE.getB()) || DE == new DEAveragedShiftedHistogramWeight(DE.getA(),DE.getB()) ){
		    log2h[s] = Num.log2(l);
		    DE.seth(l);
		    l+=0.05;
			}*/
		/*if ( DE == new DEHistogram(DE.getA(),DE.getB())){
		log2h[s] = Num.log2((model.getMax()-model.getMin())/Math.pow(2, r));
		
		DE.seth((model.getMax()-model.getMin())/Math.pow(2, r));
		DE.setM((int)Math.pow((double)2, r));
		//t+=0.0015;
		r=r+2;
		if (r ==12) 
			r=6;}
		
		else if(DE == new DEAveragedShiftedHistogram(DE.getA(),DE.getB()) || DE == new DEAveragedShiftedHistogramWeight(DE.getA(),DE.getB()) ){
			log2h[s] = Num.log2((model.getMax()-model.getMin())/(Math.pow(2, d)*32));
			
			//DE.seth((model.getMax()-model.getMin())/Math.pow(2, d));
			DE.setM((int)Math.pow((double)2, d));
			d=d+0.01;
			if (d ==10) 
				d=6;
		}
		else {
			log2h[s] = Num.log2(l);
	    DE.seth(l);
	    l+=0.05;}*/
		
		/*log2h[s] = Num.log2((model.getMax()-model.getMin())/(Math.pow(2, d)*12));
		
		//DE.seth((model.getMax()-model.getMin())/Math.pow(2, d));
		DE.setM((int)Math.pow((double)2, d)*12);
		d=d+0.001;
		if (d ==10) 
			d=6;*/
		log2h[s] = Num.log2(l);
	    DE.seth(l);
	    l+=0.001;
		
	    /*log2h[s] = Num.log2((model.getMax()-model.getMin())/(Math.pow(2, d)*32));
		
		//DE.seth((model.getMax()-model.getMin())/Math.pow(2, d));
		DE.setM((int)Math.pow((double)2, d)*32);
		d=d+0.001;
		if (d ==10) 
			d=6;*/
		
		RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statReps, data);
		integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model,  m, data, DE, numEvalPoints);
		MISE[s]=RQMCExperimentDensity.computeDensityMISE (model, m, data, DE, numEvalPoints);
		bias[s]=RQMCExperimentDensity.computeDensityBias (model,  m, data, DE, numEvalPoints);
	    log2MISE[s] = Num.log2(MISE[s]);
	    log2IV[s] = Num.log2(integVariance[s]);
	    log2ISB[s]= Num.log2(bias[s]);
	    if (displayExec) {
	    	System.out.println( " " + size[s] +" " + PrintfFormat.f(10, 5, log2IV[s]) + " " + PrintfFormat.f(10, 5, log2ISB[s]) +
		  	          " " + PrintfFormat.f(10, 5, log2MISE[s]) + "\n");
	    }
	}	
	 
  cpuTime = timer.format();	   
}
   
   /**
    * Performs an RQMC experiment with the given model for h chosen optimally in function of n , with this series of RQMC point sets.  
   * For each set in the series, computes the bias, the variance, the integrated variance (IV), the integrated square bias (ISB), the kolmogorov smirnov test (KS), the Cramer von Mises (CvM) and   the mean integrated square error (MISE), its log in base 2.
   */
   public void testMISERateOptimalU01 (MonteCarloModelDensityU01 model, int m,
			DensityEstimator DE, int numEvalPoints,  double[] MISE, double[] integVariance, double[] bias, RQMCPointSet [] theSets, String point, String densityestimator) {
	int n;
	int r = 2;
	Tally statReps = new Tally();
	
	Chrono timer = new Chrono();
	numReplicates = m;
	this.model = model;
  if (displayExec) {
  	System.out.println("\n ============================================= ");
  	System.out.println("RQMC simulation for density estimation with optimal h:  ");
  	System.out.println("Model: " + model.toString());
  	System.out.println(" Density Esimator  = " + densityestimator + "\n");
  	System.out.println(" Number of indep copies m  = " + m);
  	System.out.println(" Point sets: " + point + "\n");    
  	System.out.println("    n    log2(IV) log2(ISB)  log2(MISE)   \n");      
  }
  
 
  testMISERateVariedhnU01(model, numReplicates, DE, numEvalPoints,MISE,integVariance,bias, theSets, point); 
  
  double[] regCoefbias= slope (log2h, log2ISB, numSkipRegression);
	double alpha=regCoefbias[1];
	System.out.println("alpha"+regCoefbias[1] );

	double[] regCoeff = slope( log2n, log2h, log2IV, numSkipRegression);
	//double C=Math.exp(regCoeff[0]);
	double C=Math.pow(2, regCoeff[0]);
	//double B=Math.exp(regCoefbias[0]);
	double B=Math.pow(2,regCoefbias[0]);
	
	
	double  beta  =  -regCoeff[1] ;
	double delta =   -regCoeff[2];	
	double gamma = beta/(alpha+delta);	
	double kappa= Math.pow((C*delta/(B*alpha)),1/(alpha+delta));
	System.out.println("delta"+delta );
	System.out.println("gamma"+gamma );
	System.out.println("kappa"+kappa );
	System.out.println("B"+B );
	System.out.println("C"+C );
	System.out.println("beta"+beta );
	
  
 
	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		size[s] = n;
		double[][] data = new double[numReplicates][n];		
		log2n[s] = Num.log2(n);
		//DE.setM((int)(64*Math.pow(n,1/3)));
		//DE.seth(1/Math.pow(8, r)*Math.pow(n, -gamma));
		
		//DE.setM((int) (2000*Math.pow(n, 1/5)));
		//DE.setM((int) (64*Math.pow(n, 1/3)));
		//RQMCExperiment.simulReplicatesRQMCSaveU01 (model, theSets[s], m, statReps, statKS, statCVM, data);
		
		DE.seth(kappa*Math.pow(n, -gamma));
		
		/*if ( DE == new DEHistogram(DE.getA(),DE.getB()))
		/////// histogram  ///DE.setM((int)(64*Math.pow(n, 0.27)));
		else if (DE == new DEAveragedShiftedHistogram(DE.getA(),DE.getB()) || DE == new DEAveragedShiftedHistogramWeight(DE.getA(),DE.getB()) )
			DE.setM((int)(1024*Math.pow(n, 0.1)));
		else
			DE.seth((int)((1/1024)*Math.pow(n, -0.1)));*/
			
		//DE.setM((int)(1024*Math.pow(n,1/5)));
	     /////// histogram  ///DE.setM((int)(64*Math.pow(n, 0.27)));
		////// ASH   ///DE.setM((int)(1024*Math.pow(n, 0.20)));
		//DE.seth((1/64)*Math.pow(n, -0.20));
		/*if ( DE == new DEHistogram(DE.getA(),DE.getB())){
			DE.setM((int)(64*Math.pow(n, 0.27)));}
		else if (DE == new DEAveragedShiftedHistogram(DE.getA(),DE.getB()) || DE == new DEAveragedShiftedHistogramWeight(DE.getA(),DE.getB()) )
			DE.setM((int)(1024*Math.pow(n, 0.20)));
		else if (DE == new DEGaussianKDEBotev(DE.getA(),DE.getB()) )
			DE.setM((int)(1024*Math.pow(n, 0.20)));*/
		
		//DE.seth(kappa*Math.pow(n, -gamma));
		//DE.setM((int)(32*Math.pow(n, 0.24)));
		//DE.seth(10 / (1024.0 * Math.pow((double) n, -0.24)));
		RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], numReplicates, statReps, data);
		integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model,  numReplicates, data, DE, numEvalPoints);
		MISE[s]=RQMCExperimentDensity.computeDensityMISE (model, numReplicates, data, DE, numEvalPoints);
		bias[s]=RQMCExperimentDensity.computeDensityBias (model,  numReplicates, data, DE, numEvalPoints);
	    log2MISE[s] = Num.log2(MISE[s]);
	    log2IV[s] = Num.log2(integVariance[s]);
	    log2ISB[s]= Num.log2(bias[s]);
	    if (displayExec) {
	    	System.out.println( " " + size[s] +" " + PrintfFormat.f(10, 5, log2IV[s]) + " " + PrintfFormat.f(10, 5, log2ISB[s]) +
		  	          " " + PrintfFormat.f(10, 5, log2MISE[s]) + "\n");
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
	
	
	
	public String reportMISEVAriedhnU01 (boolean details , String densityestimator, String point) {
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for density estimation: \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + numReplicates + "\n");
		sb.append(" Density Esimator  = " + densityestimator + "\n");
		sb.append(" Point sets: " + point + "\n\n");		
		sb.append("RQMC Mean Integrated Square Error (MISE) \n");
		if (details) {
			sb.append("    n    log2(IV) log2(ISB)  log2(MISE)   \n");  
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				sb.append( " " + size[s] +" " + PrintfFormat.f(10, 5, log2IV[s]) + " " + PrintfFormat.f(10, 5, log2ISB[s]) +
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
		
	public String reportMISEU01 (boolean details, String densityestimator, String point) {
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for density estimation : \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + numReplicates + "\n");
		sb.append(" Density Esimator  = " + densityestimator + "\n");
		sb.append(" Point sets: " + this.toString() + "\n\n");
		sb.append("RQMC Mean Integrated Square Error (MISE) \n");
		if (details) {
			sb.append("    n      log2(IV)  log2(ISB) log2(MISE)   \n");  
			for (int s = 0; s < numSets; s++) { // For each cardinality n
					sb.append( " " + size[s] + " " + PrintfFormat.f(10, 5, log2IV[s]) + " " + PrintfFormat.f(10, 5, log2ISB[s]) +
				  	          " " + PrintfFormat.f(7, 2, log2MISE[s]) + "\n");
				}
			}
		
		double[] regCoeff = slope(log2n, log2MISE, numSkipRegression);
		sb.append("  Slope of log2(MISE) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");
		regCoeff = slope(log2n, log2IV, numSkipRegression);
		sb.append("  Slope of log2(IV) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");
		regCoeff = slope(log2n, log2ISB, numSkipRegression);
		sb.append("  Slope of log2(ISB) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		sb.append("  Total CPU Time = " + cpuTime + "\n");
		sb.append("-----------------------------------------------------\n");		
		return sb.toString();
	}
	
	
	/**
	 * Performs an experiment (testVarianceRate)  with the given model for each point set series in the given list,
	 * and returns a report as a string. 
	 * 
	 */
	
	
	public String TestRQMCManyPointTypes (MonteCarloModelDensityU01 model, 
			ArrayList<RQMCPointSet[]> list, int m,
			ArrayList<DensityEstimator> listDE, int numEvalPoints, 
            boolean details, String[] points) {
		StringBuffer sb = new StringBuffer("");
		numReplicates = m;	
		double[] MISE= new double[numSets];   // Will contain the IV estimates, for each n.
		double[] IntegVariance= new double[numSets];
		double[] bias= new double[numSets];
		for(int i=0; i < listDE.size(); i++) {
			int j=0;
		  for (RQMCPointSet[] ptSeries : list) {	
         	testMISERateU01 (model, m, listDE.get(i), numEvalPoints, MISE, IntegVariance,bias, ptSeries);  
  			sb.append (reportMISEU01 (details, listDE.get(i).toString(), points[j]));	
  			makePlotsMISE (numSets, m, (model.toString()).split(" ")[0], listDE.get(i).toString(), points[j]);
  			j++;
			
		  }
	    }
		
		for (RQMCPointSet[] ptSeries : list) {
			teststatistics (model,m, ptSeries);
			double[] regCoeff = slope(log2n, log2Var, numSkipRegression);
			sb.append("  Slope of log2(Var) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
			regCoeff = slope(log2n, log2KS, numSkipRegression);
			sb.append("  Slope of log2(KS) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
			regCoeff = slope(log2n, log2CVM, numSkipRegression);
			sb.append("  Slope of log2(CVM) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
			}
		return sb.toString();
	}
	public String TestRQMCManyPointTypesOptimal (MonteCarloModelDensityU01 model, 
			ArrayList<RQMCPointSet[]> list, int m,
			ArrayList<DensityEstimator> listDE, int numEvalPoints, 
            boolean details, String[] points) {
		StringBuffer sb = new StringBuffer("");
		numReplicates = m;	
		
		double[] MISE= new double[numSets];   // Will contain the IV estimates, for each n.
		double[] IntegVariance= new double[numSets];
		double[] bias= new double[numSets];
		for(int i=0; i < listDE.size(); i++) {
			int j=0;
		  for (RQMCPointSet[] ptSeries : list) {					  
         	testMISERateOptimalU01 (model, m, listDE.get(i), numEvalPoints, MISE, IntegVariance,bias, ptSeries, points[j], listDE.get(i).toString());         	
  			sb.append (reportMISEU01 (details, listDE.get(i).toString(), points[j]));	
  			makePlotsMISE (numSets, m, (model.toString()).split(" ")[0], listDE.get(i).toString(), points[j]);
  			j++;
			
		  }
	    }
		
		for (RQMCPointSet[] ptSeries : list) {
			teststatistics (model,m, ptSeries);
			double[] regCoeff = slope(log2n, log2Var, numSkipRegression);
			sb.append("  Slope of log2(Var) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
			regCoeff = slope(log2n, log2KS, numSkipRegression);
			sb.append("  Slope of log2(KS) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
			regCoeff = slope(log2n, log2CVM, numSkipRegression);
			sb.append("  Slope of log2(CVM) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
			}
		return sb.toString();
	}
	
	public String TestRQMCManyPointTypesVariedhn (MonteCarloModelDensityU01 model, 
			ArrayList<RQMCPointSet[]> list, int m,
			ArrayList<DensityEstimator> listDE, int numEvalPoints, 
            boolean details, String[] points) {
		StringBuffer sb = new StringBuffer("");
		numReplicates = m;	
		double[] MISE= new double[numSets]; //Will contain the IV estimates, for each n.
		double[] integVariance= new double[numSets];
		double[] bias= new double[numSets];
		for(int i=0; i < listDE.size(); i++) {
			int j=0;
		  for (RQMCPointSet[] ptSeries : list) {
         	testMISERateVariedhnU01 (model, m, listDE.get(i), numEvalPoints, MISE, integVariance,bias, ptSeries, points[j]);         	
			sb.append ( reportMISEVAriedhnU01 (details, listDE.get(i).toString(), points[j]));				
			makePlotsMISE (numSets, m, (model.toString()).split(" ")[0], listDE.get(i).toString(), points[j]);
			j++;
		  }
	    }
		
		for (RQMCPointSet[] ptSeries : list) {
			teststatistics (model,m, ptSeries);
			double[] regCoeff = slope(log2n, log2Var, numSkipRegression);
			sb.append("  Slope of log2(Var) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
			regCoeff = slope(log2n, log2KS, numSkipRegression);
			sb.append("  Slope of log2(KS) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
			regCoeff = slope(log2n, log2CVM, numSkipRegression);
			sb.append("  Slope of log2(CVM) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
			}
		return sb.toString();
	}
	
	 
	
	
	/** Method that produces Latex plots for the experiment for MISE
	 */
	
	public void makePlotsMISE (int numSets, int m, String descModel, String densityEstimator, String descPoints) {

		try {
			
			    chartMISE.add(log2n, log2MISE, " ", " ");
			FileWriter file = new FileWriter(descModel + "_" +   densityEstimator+ "_" + descPoints + "_MISEU01.tex");
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

/*if (DE == new DEKernelDensity(DE.getA(),DE.getB())){

log2h[s] = Num.log2(t);
DE.seth(t);
t+=0.05;	}

else
	
{
	log2h[s] = Num.log2((model.getMax()-model.getMin())/t);
	DE.seth((model.getMax()-model.getMin())/t);
	t+=0.05;	
	
}*/
/*if (DE == new DEKernelDensity(DE.getA(),DE.getB())){
log2h[s] = Num.log2(t);
DE.seth(t);
t+=0.05;}
else{
log2h[s] = Num.log2(l);
DE.seth(l);
l+=0.02;}*/
/*if ( DE == new DEHistogram(DE.getA(),DE.getB()) ){
    log2h[s] = Num.log2((DE.getB()-DE.getA())/Math.pow(2, l));
    DE.seth((DE.getB()-DE.getA())/Math.pow(2, l));
    l++;
	}*/

/*if ( DE == new DEHistogram(DE.getA(),DE.getB()) ){
    log2h[s] = Num.log2(l);
    DE.seth(l);
    l+=10;
	}*/