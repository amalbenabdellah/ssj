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

public class RQMCExperimentSeriesDensity2 extends RQMCExperimentSeries {

	ArrayList<DensityEstimator> listDE=null ;
  	public double[][] meanD = new double[3][] ;		
  	public  double[][] log2VarS = new double[3][] ;
  	double[] log2h = new double[numSets]; // log_2 of h
   /**
    * Constructor with a give series of RQMC point sets.
    *  @param theSets      the RQMC point sets
    */
   /*public RQMCExperimentSeriesDensity (RQMCPointSet[] theSets) {
	   super(theSets);
   }*/
   
   
     public RQMCExperimentSeriesDensity2 (RQMCPointSet[] theSets, ArrayList<DensityEstimator> listDE) {
  	   super(theSets);
  	   listDE =new ArrayList<DensityEstimator>();
  	   meanD = new double[3][] ;		
  	   log2VarS = new double[3][] ;
  	   this.listDE = listDE;
  	 log2h = new double[numSets];  //   log_2 of the bandwidth
     }
   public RQMCExperimentSeriesDensity2 (ArrayList<RQMCPointSet[]> theSets) {
	   super(theSets);
   }


   /**
    * Performs an RQMC experiment with the given model, with this series of RQMC point sets and a series of density estimator.  
    * For each set in the series, computes the average, the variance, its log in base 2.
    */
   
   
  /* public void testVarianceRate (MonteCarloModelBounded model, int m,
			ArrayList<DensityEstimator> listDE, int numEvalPoints, 
           double[] integVariance, RQMCPointSet [] theSets, boolean variedH) {
	   
		   for(int i=0; i<listDE.size(); i++) 		   
		   testVarianceRate (model, m, listDE.get(i), numEvalPoints, integVariance, theSets, variedH);		   
	          	  
   }*/
   
   /**
    * Performs an RQMC experiment with the given model, with this series of RQMC point sets.  
    * For each set in the series, computes the average, the variance, its log in base 2.
    */
   
   
   public void testVarianceRate (MonteCarloModelDensityKnown model, int m,
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
		//log2h[s]= -0.27*log2n[s];
		RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statReps, data);	
		
		if (VariedH==true){
		
		if(DE == new DEHistogram(DE.getA(),DE.getB())){
			DE.seth((DE.getB()-DE.getA())*Math.pow(8, r)*Math.pow(n, 0.27));
			//log2h[s] = Num.log2(DE.getB()-DE.getA()*Math.pow(8, r)*Math.pow(n, 0.27));
			
		}
		
		DE.seth(Math.pow(n, -0.27)*1/Math.pow(8, r));
		//log2h[s] = -0.27*log2n[s]+Math.log(1/Math.pow(8, r));
		//r+=3;
		}
		
		else {
		
		log2h[s] = Math.log(DE.geth());
				
		
		}
		integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model,  m, data, DE, numEvalPoints);
		//mean[s] = statReps.average();
	    log2Var[s] = Num.log2(integVariance[s]);
	    if (displayExec) {
		   System.out.println("  " + n + "     " + timer.format() + 
				      "   " + PrintfFormat.f(7, 2, log2Var[s]));
	    }
	}	
	 
  cpuTime = timer.format();	   
}
   
   public void testVarianceRate (MonteCarloModelBounded model, int m,
		   ArrayList<DensityEstimator> listDE, int numEvalPoints, 
           double[][] integVariance, RQMCPointSet [] theSets) {
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


    for(int i=0; i<listDE.size(); i++){	
		int r = 2;
	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		size[s] = n;
		double[][] data = new double[m][];
		log2n[s] = Num.log2(n);	
		/*if ( DE == new DEHistogram(DE.getA(),DE.getB()) )
		    log2h[s] = Math.log((DE.getB()-DE.getA())/Math.pow(8, l));
			 //log2h[s] = DE.geth();
		//log2h[s] = DE.geth();		
		//log2h[s]= r*Math.log(2);
		log2h[s]= Math.log(1/Math.pow(2,r));
		l++;
		r++;*/
		/*if (variedH == true){
			
			listDE.get(i).seth(Math.pow(n, -0.27)*1/Math.pow(8, r));
			r+=3;
			}*/
		//log2h[s]= Math.log(DE.geth());
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
    * 
    * 
    */
   
   
  
   
   /*public double[] regressionLogVarianceVariedh (int start) {

		
		double[] regDataX = new double[ listDE.size()*(numSets-start)];		
		double[] regDataYVar = new double[listDE.size()*(numSets-start)];

		for (int j = 0; j < listDE.size(); j++) { 			
			for (int s = 0; s < numSets-start; s++) { // For each cardinality n
				regDataX[j*(numSets-start)+s] = log2n[start+s];				
				regDataYVar[ j*(numSets-start)+s] = log2VarS[j][start+s];
			}
		}
		return  LeastSquares.calcCoefficients(regDataX, regDataYVar);
	}*/
   
   public double[] regressionLogVariance (int start) {
		double[][] regDataX = new double[ listDE.size()*(numSets-start)][2];		
		double[] regDataYMISE = new double[listDE.size()*(numSets-start)];
		
		double logh;
		for (int j = 0; j < listDE.size(); j++) { 			
				logh = Num.log2(listDE.get(j).geth());  
			for (int s = 0; s < numSets-start; s++) { // For each cardinality n
				regDataX[j*(numSets-start)+s][0] = log2n[start+s];
				regDataX[j*(numSets-start)+s][1] = logh;
				regDataYMISE[ j*(numSets-start)+s] = log2VarS[j][start+s];
			}
		}
		return  LeastSquares.calcCoefficients0(regDataX, regDataYMISE);		
	}
	public String reportVarianceFixedh (boolean details) {
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
		sb.append("  Slope of log2(var) (beta) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		sb.append("    constant term      = " + PrintfFormat.f(8, 5, Math.exp(regCoeff[0])) + "\n\n");
		sb.append("  delta = " + -regCoeff[2] + "\n");	
		sb.append("  Total CPU Time = " + cpuTime + "\n");
		sb.append("-----------------------------------------------------\n");		
		return sb.toString();
	}
	
	
	
	
	
	
/*	public double[] regressionLogMISED (int start) {
		double[][] regDataX = new double[listDE.size() * numSets][2];		
		double[] regDataYMISE = new double[listDE.size() * numSets];
		
		double logh;
		for (int j = 0; j < listDE.size(); j++) { 			
				logh = Num.log2(listDE.get(j).geth());  
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				regDataX[j  * numSets + s][0] = log2n[start+s];
				regDataX[j  * numSets + s][1] = logh;
				regDataYMISE[j * numSets + s] = log2VarS[j][start+s];
			}
		}
		return  LeastSquares.calcCoefficients0(regDataX, regDataYMISE);

	}*/
	
	
	
	
	
	/*public String reportVarianceVariedH (boolean details, double alpha) {
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
		//sb.append("  delta = " + -regCoeff[2] + "\n");		
		//sb.append("  gamma = " + (-regCoeff[1])/(alpha - regCoeff[2]) + "\n");	
		//sb.append("  nu    = " + (-alpha * regCoeff[1])/(alpha - regCoeff[2]) + "\n\n");	
		sb.append("  Total CPU Time = " + cpuTime + "\n");
		sb.append("-----------------------------------------------------\n");		
		return sb.toString();
	}*/
	
	
	
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
		double[][] integVariance= new double[listDE.size()][numSets];   // Will contain the IV estimates, for each n.
		//for(int i=0; i < listDE.size(); i++) {
		  for (RQMCPointSet[] ptSeries : list) {			
         	testVarianceRate (model, m, listDE, numEvalPoints, integVariance, ptSeries);
			sb.append (reportVarianceFixedh (details));			
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