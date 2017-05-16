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

public class RQMCExperimentSeriesDensityKnown2 extends RQMCExperimentSeriesDensity2 {


   /**
    * Constructor with a give series of RQMC point sets.
    *  @param theSets      the RQMC point sets
    */
	double[] log2h;  //   log_2 of the bandwidth
	double[] log2MISE;  
	 ArrayList<DensityEstimator> listDE=null ;
	public double[][] log2MISES = new double[3][] ;		
	public  double[][] log2VarS = new double[3][] ;
   public RQMCExperimentSeriesDensityKnown2 (RQMCPointSet[] theSets, ArrayList<DensityEstimator> listDE) {
	   super(theSets, listDE);
	   listDE =new ArrayList<DensityEstimator>();
	   log2MISES = new double[3][] ;		
	   log2VarS = new double[3][] ;
	   log2MISE = new double[numSets];
	   this.listDE = listDE;
	   log2h = new double[numSets];  //   log_2 of the bandwidth
   }
   


   /**
    * Performs an RQMC experiment with the given model, with this series of RQMC point sets and a series of density estimator.  
    * For each set in the series, computes the average, the variance, its log in base 2.
    */
   
   
   /*public void testMISERate (MonteCarloModelDensityKnown model, int m,
			ArrayList<DensityEstimator> listDE, int numEvalPoints, 
           double[] MISE, double[] integVariance, RQMCPointSet [] theSets) {
	   
		   for(int i=0; i<listDE.size(); i++) 		   
		  // testMISERate (model, m, listDE.get(i), numEvalPoints, MISE, theSets);	
			   testMISERate (model, m, listDE.get(i), numEvalPoints,  MISE, integVariance,  theSets);
	          	  
   }*/
   
   
   
   public void testMISERate (MonteCarloModelDensityKnown model, int m,
			DensityEstimator DE, int numEvalPoints,  double[] MISE, double[] integVariance,  RQMCPointSet [] theSets, boolean VariedH) {
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

  int r = 2;
 
	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		size[s] = n;
		double[][] data = new double[m][];
		log2n[s] = Num.log2(n);
		RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statReps, data);	
		
		if (VariedH==true){
		
		if(DE == new DEHistogram(DE.getA(),DE.getB())){
			DE.seth((DE.getB()-DE.getA())*Math.pow(8, r)*Math.pow(n, 0.27));
		}
		
		DE.seth(Math.pow(n, -0.27)*1/Math.pow(8, r));
		}
		
		else {
		
		log2h[s] = Math.log(DE.geth());
				
		
		}
		integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model,  m, data, DE, numEvalPoints);
		MISE[s]=RQMCExperimentDensity.computeDensityMISE (model, m, data, DE, numEvalPoints);
		//mean[s] = statReps.average();
	    log2MISE[s] = Num.log2(MISE[s]);
	    log2Var[s] = Num.log2(integVariance[s]);
	    if (displayExec) {
		   System.out.println("  " + n + "     " + timer.format() + 
				      "   " + PrintfFormat.f(7, 2, log2MISE[s]));
	    }
	}	
	 
  cpuTime = timer.format();	   
}
   
   
   public void testMISERateD (MonteCarloModelDensityKnown model, int m,
		   ArrayList<DensityEstimator> listDE, int numEvalPoints,  double[][] MISE, double[][] integVariance,  RQMCPointSet [] theSets, boolean VariedH) {
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

     log2MISES= new double[listDE.size()][numSets];
     log2VarS = new double[listDE.size()][numSets];
     
     int r = 2;
  
	for(int i=0; i<listDE.size(); i++){	
		
 
	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
	
		size[s] = n;
		double[][] data = new double[m][];
		log2n[s] = Num.log2(n);
		//log2h[s]= -0.27*log2n[s];
		RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statReps, data);	
		
		/*if(listDE.get(i) == new DEHistogram(listDE.get(i).getA(),listDE.get(i).getB())){
			listDE.get(i).seth((listDE.get(i).getB()-listDE.get(i).getA())*Math.pow(8, r)*Math.pow(n, 0.27));
			log2h[s] = Num.log2(listDE.get(i).getB()-listDE.get(i).getA()*Math.pow(8, r)*Math.pow(n, 0.27));
			
		}*/
		//log2h[s] = -0.27*log2n[s]+Math.log(1/Math.pow(8, r));
		if (VariedH == true){
		
		listDE.get(i).seth(Math.pow(n, -0.27)*1/Math.pow(8, r));
		
		}
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
	r=+1;
  }
	 
  cpuTime = timer.format();	   
}
   
   
 /*  public void testMISERateVariedH (MonteCarloModelDensityKnown model, int m,
		   ArrayList<DensityEstimator> listDE, int numEvalPoints,  double[][] MISE, double[][] integVariance,  RQMCPointSet [] theSets) {
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

  int r = 2;
  
  for(int i=0; i<listDE.size(); i++){	
 
	for (int s = 0; s < numSets; s++) { // For each cardinality n
		
		n = theSets[s].getNumPoints();
	
		size[s] = n;
		double[][] data = new double[m][];
		log2n[s] = Num.log2(n);
		//log2h[s]= -0.27*log2n[s];
		RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statReps, data);	
		
		if(listDE.get(i) == new DEHistogram(listDE.get(i).getA(),listDE.get(i).getB())){
			listDE.get(i).seth((listDE.get(i).getB()-listDE.get(i).getA())*Math.pow(8, r)*Math.pow(n, 0.27));
			log2h[s] = Num.log2(listDE.get(i).getB()-listDE.get(i).getA()*Math.pow(8, r)*Math.pow(n, 0.27));
			
		}
	
		
		listDE.get(i).seth(Math.pow(n, -0.27)*1/Math.pow(8, r));
		//log2h[s] = -0.27*log2n[s]+Math.log(1/Math.pow(8, r));
		r++;
		integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model, m, data, listDE.get(i), numEvalPoints);
		MISE[s]=RQMCExperimentDensity.computeDensityMISE (model, m, data, listDE.get(i), numEvalPoints);
		//mean[s] = statReps.average();
	    log2MISE[s] = Num.log2(MISE[s]);
	    log2Var[s] = Num.log2(integVariance[s]);
	    if (displayExec) {
		   System.out.println("  " + n + "     " + timer.format() + 
				      "   " + PrintfFormat.f(7, 2, log2MISE[s]));
	    }
	}	
  }
	 
  cpuTime = timer.format();	   
}*/
   /**
    * Performs an RQMC experiment with the given model, with this series of RQMC point sets.  
    * For each set in the series, computes the average, the variance, its log in base 2.
    */
   
   
  /* public void testMISERate (MonteCarloModelDensityKnown model, int m,
			DensityEstimator DE, int numEvalPoints,  double[] MISE, double[] integVariance,  RQMCPointSet [] theSets) {
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

   int r = 2;
  
	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		size[s] = n;
		double[][] data = new double[m][];
		log2n[s] = Num.log2(n);
		//log2h[s]= -0.27*log2n[s];
		RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statReps, data);	
		
		if(DE == new DEHistogram(DE.getA(),DE.getB())){
			DE.seth((DE.getB()-DE.getA())*Math.pow(8, r)*Math.pow(n, 0.27));
			log2h[s] = Num.log2(DE.getB()-DE.getA()*Math.pow(8, r)*Math.pow(n, 0.27));
			
		}
	
		
		DE.seth(Math.pow(n, -0.27)*1/Math.pow(8, r));
		log2h[s] = -0.27*log2n[s]+Math.log(1/Math.pow(8, r));
		r++;
		integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model, n, m, data, DE, numEvalPoints);
		MISE[s]=RQMCExperimentDensityKnown.computeDensityMISE (model, n, m, data, DE, numEvalPoints);
		//mean[s] = statReps.average();
	    log2MISE[s] = Num.log2(MISE[s]);
	    log2Var[s] = Num.log2(integVariance[s]);
	    if (displayExec) {
		   System.out.println("  " + n + "     " + timer.format() + 
				      "   " + PrintfFormat.f(7, 2, log2MISE[s]));
	    }
	}	
	 
   cpuTime = timer.format();	   
}*/



   /**
    * Returns the vector of log_2(n).
    */
   public double[] getLog2n() {
      return log2n;
   }

   
   /*public double[] regressionLogMISEVariedh (int numSkip) {
		double[] x2 = new double[numSets-numSkip], y2 = new double[numSets-numSkip];
		
		for (int j = 0; j <listDE.size(); ++j){
		for (int i = 0; i < numSets-numSkip; ++i) {
			x2[i] = log2n[i+numSkip];
			y2[i] = log2MISES[j][i+numSkip];
		}}
		return LeastSquares.calcCoefficients(x2, y2, 1);
	}*/
   
   
   
   /*public double[] regressionLogMISEVariedh (int start) {
		double[] regDataX = new double[listDE.size() * numSets];
		double[] regDataY = new double[listDE.size()  *numSets];
		double[] regDataYMISE = new double[listDE.size() * numSets];
		double[] coef;
		double logh;
		for (int j = 0; j < listDE.size(); j++) { 			
				//logh = Num.log2(listDE.get(j).geth());  
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				regDataX[j  * numSets + s] = log2n[start+s];
				regDataYMISE[j * numSets + s] = log2MISES[j][start+s];
			}
		}
		return  LeastSquares.calcCoefficients(regDataX, regDataYMISE, 1);

	}*/
   
   /*public double[] regressionLogMISED (int start) {
		double[][] regDataX = new double[listDE.size() * numSets-start][2];		
		double[] regDataYMISE = new double[listDE.size() * numSets-start];
		
		double logh;
		for (int j = 0; j < listDE.size(); j++) { 			
				logh = Num.log2(listDE.get(j).geth());  
			for (int s = 0; s < numSets-start; s++) { // For each cardinality n
				regDataX[j  * numSets + s][0] = log2n[start+s];
				regDataX[j  * numSets + s][1] = logh;
				regDataYMISE[j * numSets + s] = log2MISES[j][start+s];
			}
		}
		return  LeastSquares.calcCoefficients0(regDataX, regDataYMISE);

	}*/
   
   /*public double[] regressionLogMISED (int start) {
		double[][] regDataX = new double[ numSets-start][2];		
		double[] regDataYMISE = new double[numSets-start];
		
		double logh;
		for (int j = 0; j < listDE.size(); j++) { 			
				logh = Num.log2(listDE.get(j).geth());  
			for (int s = 0; s < numSets-start; s++) { // For each cardinality n
				regDataX[j  * numSets + s][0] = log2n[start+s];
				regDataX[j  * numSets + s][1] = logh;
				regDataYMISE[j * numSets + s] = log2MISES[j][start+s];
			}
		}
		return  LeastSquares.calcCoefficients0(regDataX, regDataYMISE);

	}
   */
   
   /*public double[] regressionLogMISED (int start) {
		double[][] regDataX = new double[ numSets-start][2];		
		double[] regDataYMISE = new double[numSets-start];
		
		double logh;
		for (int j = 0; j < listDE.size(); j++) { 			
				logh = Num.log2(listDE.get(j).geth());  
			for (int s = 0; s < numSets-start; s++) { // For each cardinality n
				regDataX[s][0] = log2n[start+s];
				regDataX[ s][1] = logh;
				regDataYMISE[ s] = log2MISES[j][start+s];
			}
		}
		return  LeastSquares.calcCoefficients0(regDataX, regDataYMISE);
		

	}*/
   
   public double[] regressionLogMISED (int start) {
		double[][] regDataX = new double[ listDE.size()*(numSets-start)][2];		
		double[] regDataYMISE = new double[listDE.size()*(numSets-start)];
		
		double logh;
		for (int j = 0; j < listDE.size(); j++) { 			
				logh = Num.log2(listDE.get(j).geth());  
			for (int s = 0; s < numSets-start; s++) { // For each cardinality n
				regDataX[j*(numSets-start)+s][0] = log2n[start+s];
				regDataX[j*(numSets-start)+s][1] = logh;
				regDataYMISE[ j*(numSets-start)+s] = log2MISES[j][start+s];
			}
		}
		return  LeastSquares.calcCoefficients0(regDataX, regDataYMISE);
		

	}
   
   /*public double[] regressionLogMISEVariedh (int start) {

		
		double[][] regDataX = new double[ listDE.size()*(numSets-start)][1];		
		double[] regDataYMISE = new double[listDE.size()*(numSets-start)];

		for (int j = 0; j < listDE.size(); j++) { 			
			for (int s = 0; s < numSets-start; s++) { // For each cardinality n
				regDataX[j*(numSets-start)+s][0] = log2n[start+s];				
				regDataYMISE[ j*(numSets-start)+s] = log2MISES[j][start+s];
			}
		}
		return  LeastSquares.calcCoefficients(regDataX, regDataYMISE);
	}*/
   
   public double[] regressionLogMISEVariedh (int start) {

		
		double[] regDataX = new double[ listDE.size()*(numSets-start)];		
		double[] regDataYMISE = new double[listDE.size()*(numSets-start)];

		for (int j = 0; j < listDE.size(); j++) { 			
			for (int s = 0; s < numSets-start; s++) { // For each cardinality n
				regDataX[j*(numSets-start)+s] = log2n[start+s];				
				regDataYMISE[ j*(numSets-start)+s] = log2MISES[j][start+s];
			}
		}
		return  LeastSquares.calcCoefficients(regDataX, regDataYMISE,0);
	}
  
   
   /**
    * Produces and returns a report on the last experiment.
    * @param numSkip  The first numSkip values of n are skipped for the regression
    * @param details  If true, gives values (mean, log variance,...) for each n.
    * @return  Report as a string.
    */
	public String reportMISEVariedH (boolean details, ArrayList<DensityEstimator> listDE) {
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for density estimation, with known density: \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + numReplicates + "\n");
		sb.append(" Point sets: " + this.toString() + "\n\n");
		sb.append("RQMC Mean Integrated Square Error (MISE) \n");
		if (details) {
								
			for (int i = 0;  i< listDE.size(); i++) { 
				sb.append("    n  log2(Var)  log2(MISE) \n");
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				sb.append(" " + size[s] + " " + PrintfFormat.f(10, 5, log2VarS[i][s]) +
				          " " + PrintfFormat.f(7, 2, log2MISES[i][s]) + "\n");
			}
			}
		}
		double[] regCoeff = regressionLogMISEVariedh (numSkipRegression);
		sb.append("  Slope of log2(MISE) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");
		sb.append("  Total CPU Time = " + cpuTime + "\n");
		sb.append("-----------------------------------------------------\n");		
		return sb.toString();
	}
	
	
	/*public double[] regressionLogMISE (int numSkip, ArrayList<DensityEstimator> listDE) {
		double[][] x2 = new double[numSets-numSkip][2];
		double [] y2 = new double[numSets-numSkip];
		for(int j=0; j<listDE.size(); j++){
		for (int i = 0; i < numSets-numSkip; ++i) {
			x2[i][0] = log2n[i+numSkip];
			x2[i][1] = Num.log2(listDE.get(j).geth());			
			y2[i] = log2MISE[i+numSkip];
		}}
		return LeastSquares.calcCoefficients0(x2, y2);
	}*/
	
	public String reportMISEFixedH (boolean details, double alpha, ArrayList<DensityEstimator> listDE) {
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for density estimation, with known density: \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + numReplicates + "\n");
		sb.append(" Point sets: " + this.toString() + "\n\n");
		sb.append("RQMC Mean Integrated Square Error (MISE) \n");
		if (details) {
			for (int i = 0;  i< listDE.size(); i++) { 
				sb.append("    n  log2(Var)  log2(MISE) \n");
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				sb.append(" " + size[s] + " " + PrintfFormat.f(10, 5, log2VarS[i][s]) +
				          " " + PrintfFormat.f(7, 2, log2MISES[i][s]) + "\n");
			}
			}
		}
		double[] regCoeff = regressionLogMISED(numSkipRegression);
		//sb.append("  Slope of log2(var) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		//sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");
		sb.append("  Constant   for MISE  = " + regCoeff[0] + "\n");
		sb.append("  beta for MISE = " + -regCoeff[1] + "\n");
		double delta = alpha - regCoeff[2];
		sb.append("  delta for MISE = " + delta + "\n");			
		sb.append("  gamma = " + (-regCoeff[1])/(alpha - regCoeff[2]) + "\n");	
		sb.append("  nu    = " + (-alpha * regCoeff[1])/(alpha - regCoeff[2]) + "\n\n");
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
		double[][] MISE= new double[listDE.size()][numSets]; //Will contain the IV estimates, for each n.
		double[][] integVariance= new double[listDE.size()][numSets];
		//for(int i=0; i < listDE.size(); i++) {
		  for (RQMCPointSet[] ptSeries : list) {			
         	testMISERateD (model, m, listDE, numEvalPoints, MISE, integVariance, ptSeries, true);
         	for (DensityEstimator lDE : listDE) {
         	if ( lDE.equals(new DEHistogram(model.getMin(),model.getMax())))
  			  sb.append (reportMISEFixedH (details,2,listDE));	
           	else          		
           	  sb.append (reportMISEFixedH (details,4, listDE));
			
		  }
		  }
		return sb.toString();
	}
	
	public String TestRQMCManyPointTypes (MonteCarloModelDensityKnown model, 
			ArrayList<RQMCPointSet[]> list, int m,
			ArrayList<DensityEstimator> listDE, int numEvalPoints, 
            boolean details, boolean VariedH) {
		StringBuffer sb = new StringBuffer("");
		numReplicates = m;	
		double[][] MISE= new double[listDE.size()][numSets]; //Will contain the IV estimates, for each n.
		double[][] integVariance= new double[listDE.size()][numSets];
		  for (RQMCPointSet[] ptSeries : list) {			
         	testMISERateD (model, m, listDE, numEvalPoints, MISE, integVariance, ptSeries, VariedH);   
         	if (VariedH == true)
			sb.append ( reportMISEVariedH (details, listDE));
         	else {         		
         		if ( listDE.equals(new DEHistogram(model.getMin(),model.getMax())))
        			  sb.append (reportMISEFixedH (details,2,listDE));	
                 	else          		
                 	  sb.append (reportMISEFixedH (details,4, listDE));
         		         	
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