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

public class RQMCExperimentSeriesDensityKnown extends RQMCExperimentSeriesDensity {


   /**
    * Constructor with a give series of RQMC point sets.
    *  @param theSets      the RQMC point sets
    */
	XYLineChart chart = new XYLineChart ();
	
double[] log2MISE = new double[numSets]; // log_2 of MISE
public double[][] meanD ;
public double[][] log2MISES;

	
	double[] log2bias = new double[numSets]; // log_2 of h
   public RQMCExperimentSeriesDensityKnown (RQMCPointSet[] theSets) {
	   super(theSets);
	   log2MISE = new double[numSets]; // log_2 of the MISE
	   log2bias = new double[numSets];  //   log_2 of the bias
	   
   }
   public RQMCExperimentSeriesDensityKnown (ArrayList<RQMCPointSet[]> theSets) {
	   super(theSets);
	   log2MISE = new double[numSets]; // log_2 of the MISE
	   log2bias = new double[numSets];  //   log_2 of the bias
	   
   }

   
   
   /**
    * Performs an RQMC experiment with the given model, with this series of RQMC point sets and a series of density estimator.  
    * For each set in the series, computes the average, the variance, its log in base 2.
    */
   
   public void testMISERateD (MonteCarloModelDensityKnown model, int m,
		   ArrayList<DensityEstimator> listDE, int numEvalPoints,  double[][] MISE, double[][] integVariance,  RQMCPointSet [] theSets) {
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
	System.out.println("    n     CPU time         mean       log2(var) ");	    	
  }

     log2MISES= new double[listDE.size()][numSets];
     log2VarS = new double[listDE.size()][numSets];   
	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		for(int i=0; i<listDE.size(); i++){	
		size[s] = n;
		double[][] data = new double[m][];
		
		log2n[s] = Num.log2(n);		
		RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statReps, data);	
		//KS[s] = statKS.average();
		//CVM[s] = statCVM.average();
		
		integVariance[i][s]=RQMCExperimentDensity.computeDensityVariance (model, m, data, listDE.get(i), numEvalPoints);
		MISE[i][s]=RQMCExperimentDensity.computeDensityMISE (model, m, data, listDE.get(i), numEvalPoints);
		//mean[s] = statReps.average();
	    log2MISES[i][s] = Num.log2(MISE[i][s]);
	    log2VarS[i][s] = Num.log2(integVariance[i][s]);
	    if (displayExec) {
		   System.out.println("  " + n + "     " + timer.format() + PrintfFormat.f(7, 2,log2VarS[i][s])+
				      "   " + PrintfFormat.f(7, 2, log2MISES[i][s]));
	    }
	}	
  }
	 
  cpuTime = timer.format();	   
}
   /**
    * Performs an RQMC experiment with the given model, with this series of RQMC point sets.  
    * For each set in the series, computes the bias, the variance, the MISE, its log in base 2.
    */
   
   
   public void testMISERate (MonteCarloModelDensityKnown model, int m,
			DensityEstimator DE, int numEvalPoints,  double[] MISE, double[] integVariance, double[] bias, RQMCPointSet [] theSets) {
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
	System.out.println("    n     CPU time         mean       log2(var) ");	    	
   }

   int r = 2;
  
	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		size[s] = n;
		double[][] data = new double[m][];
		
		log2n[s] = Num.log2(n);
		//log2h[s]= -0.27*log2n[s];
		RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statReps,data);
		integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model,  m, data, DE, numEvalPoints);
		MISE[s]=RQMCExperimentDensity.computeDensityMISE (model, m, data, DE, numEvalPoints);
		bias[s]=RQMCExperimentDensity.computeDensityBias (model,  m, data, DE, numEvalPoints);
		//bias[s]=RQMCExperimentDensity.computeDensityMISE (model,  m, data, DE, numEvalPoints)-variance[s];
		//KS[s] = statKS.average();
		//CVM[s] = statCVM.average();
		//mean[s] = statReps.average();
	    log2MISE[s] = Num.log2(MISE[s]);
	    log2IV[s] = Num.log2(integVariance[s]);
	    log2bias[s]= Num.log2(bias[s]);
	    if (displayExec) {
		   System.out.println("  " + n + "     " + timer.format() + 
				      "   " + PrintfFormat.f(7, 2, log2MISE[s]));
	    }
	}	
	 
   cpuTime = timer.format();	   
}


   public void testMISERateVariedhn (MonteCarloModelDensityKnown model, int m,
			DensityEstimator DE, int numEvalPoints,  double[] MISE, double[] integVariance, double[] bias, RQMCPointSet [] theSets) {
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
	System.out.println("    n     CPU time         mean       log2(var) ");	    	
  }

  int r = 2;
  double  l=1;
		  // t=-5;
 
	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		size[s] = n;
		double[][] data = new double[m][];
	
		log2n[s] = Num.log2(n);
		//log2h[s]= -0.27*log2n[s];
		RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statReps, data);	
		
		
		
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
			/*log2h[s] =t*Math.log(2);
		    DE.seth(t* Math.log(2));	
		    t++;*/
		}
		
		
		
		integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model,  m, data, DE, numEvalPoints);
		MISE[s]=RQMCExperimentDensity.computeDensityMISE (model, m, data, DE, numEvalPoints);
		bias[s]=RQMCExperimentDensity.computeDensityBias (model,  m, data, DE, numEvalPoints);
		//bias[s]=RQMCExperimentDensity.computeDensityMISE (model,  m, data, DE, numEvalPoints)-variance[s];
		//bias[s] = statReps.average();
		//mean[s] = statReps.average();
		//KS[s] = statKS.average();
		//CVM[s] = statCVM.average();
	    log2MISE[s] = Num.log2(MISE[s]);
	    log2IV[s] = Num.log2(integVariance[s]);
	    log2bias[s]= Num.log2(bias[s]);
	    if (displayExec) {
		   System.out.println("  " + n + "     " + timer.format() + 
				      "   " + PrintfFormat.f(7, 2, log2MISE[s]));
	    }
	}	
	 
  cpuTime = timer.format();	   
}
   
   public void testMISERateOptimal (MonteCarloModelDensityKnown model, int m,
			DensityEstimator DE, int numEvalPoints,  double[] MISE, double[] integVariance, double[] bias, RQMCPointSet [] theSets) {
	int n;
	int r=4;
	Tally statReps = new Tally();
	Tally statKS = new Tally();
	Tally statCVM = new Tally();
	Chrono timer = new Chrono();
	numReplicates = m;
	this.model = model;
  if (displayExec) {
  	System.out.println("\n ============================================= ");
  	System.out.println("RQMC simulation for density estimation, for known density:  ");
  	System.out.println("Model: " + model.toString());
  	System.out.println(" Number of indep copies m  = " + m);
  	System.out.println(" Point sets: " + theSets[0].toString() + "\n");
	System.out.println("    n     CPU time         mean       log2(var) ");	    	
  }
  
  /*if (DE==new DEHistogram(DE.getA(),DE.getB()))
	  alpha=2;*/
  
 // testVarianceRateVariedhn(model, m, DE, numEvalPoints,integVariance, theSets); 
  testMISERateVariedhn(model, m, DE, numEvalPoints,MISE,integVariance,bias, theSets); 
  double[] regCoefbias= regressionLogBias (numSkipRegression);
	double alpha=regCoefbias[1];
	System.out.println("alpha"+regCoefbias[1] );
  
//double alpha=4;
	double[] regCoeff = regressionLogVarianceVariedhn(numSkipRegression);
	double C=Math.exp(regCoeff[0]);
	//double B=C;
	double B=Math.exp(regCoefbias[0]);
	
	
	double  beta  =  -regCoeff[1] ;
	double delta =   -regCoeff[2];	
	double gamma = beta/(alpha+delta);	
	double kappa= Math.pow((C*delta/B*alpha),1/(alpha+delta));
	System.out.println("gamma"+gamma );
  /*testMISERateVariedhn (model, m, DE, numEvalPoints, MISE, integVariance,bias, theSets); 
  double[] regCoefbias= regressionLogBias (numSkipRegression);
	double alpha=regCoefbias[1];
	double[] regCoeff = regressionLogMISEVariedhn(numSkipRegression);
	double beta =-regCoeff[1] ;
	double delta = alpha - regCoeff[2];
			
	double gamma = -regCoeff[1]/(alpha+delta);	
	double B=Math.exp(regCoefbias[0]);
	double C=Math.exp(regCoeff[0]-regCoefbias[0]);
	double kappa= Math.pow((C*delta/B*alpha),1/(alpha+delta));*/
	//double nu    = (-alpha * regCoeff[1])/(alpha - regCoeff[2]);
	
  
 
	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		size[s] = n;
		double[][] data = new double[m][];
		/*double[][] KS = new double[m][];
		double[][] CVM = new double[m][];*/
		log2n[s] = Num.log2(n);
		//log2h[s]= -0.27*log2n[s];
		
		/*if(DE == new DEHistogram(DE.getA(),DE.getB())){
			DE.seth((DE.getB()-DE.getA())*Math.pow(8, r)*Math.pow(n, 0.27));
		}*/
		
		//DE.seth(Math.pow(n, -0.27)*1/Math.pow(8, r));
		//DE.seth(kappa*Math.pow(n, -gamma));
		DE.seth(1/Math.pow(8, r)*Math.pow(n, -gamma));
		RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statReps, data);	
		integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model,  m, data, DE, numEvalPoints);
		MISE[s]=RQMCExperimentDensity.computeDensityMISE (model, m, data, DE, numEvalPoints);
		bias[s]=RQMCExperimentDensity.computeDensityBias (model,  m, data, DE, numEvalPoints);
		//KS[s] = statKS.average();
		//CVM[s] = statCVM.average();
		System.out.println("MI"+MISE[s] );
		//mean[s] = statReps.average();
	    log2MISE[s] = Num.log2(MISE[s]);
	    log2IV[s] = Num.log2(integVariance[s]);
	    //log2bias[s]= Num.log2(bias[s]);
	    if (displayExec) {
		   System.out.println("  " + n + "     " + timer.format() + 
				      "   " + PrintfFormat.f(7, 2, log2MISE[s]));
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

   public double[] regressionLogBias (int numSkip) {
		double[] x2 = new double[numSets-numSkip];
		double [] y2 = new double[numSets-numSkip];
		for (int i = 0; i < numSets-numSkip; ++i) {
			x2[i] = log2h[i+numSkip];	
			y2[i] = log2bias[i+numSkip];
		}
		return LeastSquares.calcCoefficients(x2, y2,1);
	}
   
   public double[] regressionLogMISEVariedhn (int numSkip) {
		double[][] x2 = new double[numSets-numSkip][2];
		double [] y2 = new double[numSets-numSkip];
		for (int i = 0; i < numSets-numSkip; ++i) {
			x2[i][0] = log2n[i+numSkip];
			x2[i][1] = log2h[i+numSkip];			
			y2[i] = log2MISE[i+numSkip];
		}
		return LeastSquares.calcCoefficients0(x2, y2);
	}
   
   /**
    * Produces and returns a report on the last experiment.
    * @param numSkip  The first numSkip values of n are skipped for the regression
    * @param details  If true, gives values (mean, log variance,...) for each n.
    * @return  Report as a string.
    */
	public String reportMISEVariedhn (boolean details) {
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for density estimation, with known density: \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + numReplicates + "\n");
		sb.append(" Point sets: " + this.toString() + "\n\n");
		sb.append("RQMC Mean Integrated Square Error (MISE) \n");
		if (details) {
			
			sb.append("    n    log2(IV)  log2(MISE) \n");
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				sb.append(" " + size[s] + " " + PrintfFormat.f(10, 5, log2IV[s]) +
				          " " + PrintfFormat.f(7, 2, log2MISE[s]) + "\n");
			}
		}
		double[] regCoeff = regressionLogMISEVariedhn (numSkipRegression);
		sb.append("  Slope of log2(MISE) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");
		sb.append("  Total CPU Time = " + cpuTime + "\n");
		sb.append("-----------------------------------------------------\n");		
		return sb.toString();
	}
	
	
	/*public double[] regressionLogMISED (int start, ArrayList<DensityEstimator> listDE) {
		double[][] regDataX = new double[ numSets-start][2];		
		double[] regDataYMISE = new double[numSets-start];
		
		double logh;
		for (int j = 0; j < listDE.size(); j++) { 	
				logh = Num.log2(listDE.get(j).geth());  
			for (int s = 0; s < numSets-start; s++) { // For each cardinality n
				regDataX[s][0] = log2n[start+s];
				regDataX[ s][1] = logh;
				regDataYMISE[ s] = log2MISE[start+s];
			}
		}
		return  LeastSquares.calcCoefficients0(regDataX, regDataYMISE);
		

	}*/
	
	
	
	public double[] regressionLogMISE (int numSkip) {
		double[] x2 = new double[numSets-numSkip];
		double [] y2 = new double[numSets-numSkip];
		for (int i = 0; i < numSets-numSkip; ++i) {
			x2[i] = log2n[i+numSkip];			
			y2[i] = log2MISE[i+numSkip];
		}
		return LeastSquares.calcCoefficients(x2, y2, 1);
	}
	
	public String reportMISEVAriedhn (boolean details) {
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for density estimation, with known density: \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + numReplicates + "\n");
		sb.append(" Point sets: " + this.toString() + "\n\n");
		sb.append("RQMC Mean Integrated Square Error (MISE) \n");
		if (details) {
			sb.append("    n    log2(IV)  log2(MISE) \n");
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				sb.append(" " + size[s] + " " + PrintfFormat.f(10, 5, log2IV[s]) +
				          " " + PrintfFormat.f(7, 2, log2MISE[s]) + "\n");
			}
		}
		double[] regCoefbias= regressionLogBias (numSkipRegression);
		double alpha=regCoefbias[1];
		double[] regCoeff = regressionLogMISEVariedhn(numSkipRegression);
		//sb.append("  Slope of log2(var) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		//sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");
		//sb.append("  C   for MISE  = " + Math.exp(regCoeff[0]) + "\n");
		sb.append("  beta for MISE = " + -regCoeff[1] + "\n");
		double delta = alpha - regCoeff[2];
		sb.append("  delta for MISE = " + delta + "\n");			
		sb.append("  gamma = " + (-regCoeff[1])/(alpha - regCoeff[2]) + "\n");	
		sb.append("  nu    = " + (-alpha * regCoeff[1])/(alpha - regCoeff[2]) + "\n\n");
		sb.append("  Total CPU Time = " + cpuTime + "\n");
		sb.append("-----------------------------------------------------\n");		
		return sb.toString();
	}
	public String reportMISE (boolean details) {
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for density estimation, with known density: \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + numReplicates + "\n");
		sb.append(" Point sets: " + this.toString() + "\n\n");
		sb.append("RQMC Mean Integrated Square Error (MISE) \n");
		if (details) {
			sb.append("    n    log2(IV)  log2(MISE) \n");
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				sb.append(" " + size[s] + " " + PrintfFormat.f(10, 5, log2IV[s]) +
				          " " + PrintfFormat.f(7, 2, log2MISE[s]) + "\n");
			}
		}
		double[] regCoeff = regressionLogMISE(numSkipRegression);
		sb.append("  Slope of log2(MISE) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
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
  			sb.append (reportMISE (details));	
  			makePlotsMISE (numSets, m, (model.toString()).split(" ")[0], " ");
			
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
  			sb.append (reportMISE (details));	
  			makePlotsMISE (numSets, m, (model.toString()).split(" ")[0], " ");
			
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
			sb.append ( reportMISEVariedhn (details));	
			makePlotsMISE (numSets, m, (model.toString()).split(" ")[0], " ");
		  }
	    }
		return sb.toString();
	}
	
	
	/**
	 * Performs an experiment (testVarianceRate) for each point set series in the given list,
	 * and returns a report as a string. 
	 * 
	 */
	
	public void makePlotsMISE (int numSets, int m, String descModel, String descPoints) {
		// makeGraph();
		try {
			// String title = descModel + "; " + descPoints;
			// double [][] bidon = new double[0][0];
			// XYLineChart chart = new XYLineChart (title, "lg n", "lg MISE", bidon);
			// chart.init (title, "lg n", "lg MISE");
		    
			//for (int j = 3; j < numStats; j++)
				// chart.add(log2n, log2StatsMISE[j], statNames[j], " ");
			    chart.add(log2n, log2MISE, " ", " ");
			FileWriter file = new FileWriter(descModel + "_" + descPoints + "_MISE.tex");
			file.write(chart.toLatex(12, 8));
			file.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void fitPrintRegressionMISE (int order,ArrayList<DensityEstimator> listDE, boolean useBandwidth, 
			  int start, int number,
			  double range, double alpha, String method, StringBuffer sb) {
		double[][] regDataX = new double[listDE.size() * number][order];
		double[] regDataY = new double[listDE.size()  * number];
		double[] regDataYMISE = new double[listDE.size()  * number];
		double[] coef;
		double logh;
		for (int j = 0; j < listDE.size() ; j++) { // For each m.
				logh = Num.log2(listDE.get(j).geth());  // For KDE.	
			for (int s = 0; s < number; s++) { // For each cardinality n
				regDataX[j * number + s][0] = log2n[start+s];
				regDataX[j  * number + s][1] = logh;
				if (order > 2) regDataX[j * number + s][2] = logh * log2n[start+s];
				regDataYMISE[j * number + s] = log2MISES[j][start+s];
			}
		}
		coef = LeastSquares.calcCoefficients0(regDataX, regDataY);
		// System.out.println(" coef computed 1.\n");
		sb.append("  Regression coefficients for " + method + ".\n");
		sb.append("  C     = " + Math.exp(coef[0]) + "\n");
		sb.append("  beta  = " + -coef[1] + "\n");
		sb.append("  delta = " + -coef[2] + "\n");
		if (order > 2) sb.append("  inter = " + -coef[3] + "\n");
		sb.append("  gamma = " + (-coef[1])/(alpha - coef[2]) + "\n");	
		sb.append("  nu    = " + (-alpha * coef[1])/(alpha - coef[2]) + "\n\n");	

			coef = LeastSquares.calcCoefficients0(regDataX, regDataYMISE);
			sb.append("  C for MISE     = " + Math.exp(coef[0]) + "\n");
			sb.append("  beta for MISE  = " + -coef[1] + "\n");
			sb.append("  delta for MISE = " + -coef[2] + "\n\n");
		}
	
	
	public String toString () {
		return theSets[0].toString();
	}
}