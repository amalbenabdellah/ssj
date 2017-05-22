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

import umontreal.ssj.hups.*;
import umontreal.ssj.stat.Tally;
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

public class RQMCExperimentSeriesDensityF extends RQMCExperimentSeries {


   /**
    * Constructor with a give series of RQMC point sets.
    *  @param theSets      the RQMC point sets
    */
   public RQMCExperimentSeriesDensityF (RQMCPointSet[] theSets) {
	   super(theSets);
   }
   public RQMCExperimentSeriesDensityF (ArrayList<RQMCPointSet[]> theSets) {
	   super(theSets);
   }


   /**
    * Performs an RQMC experiment with the given model, with this series of RQMC point sets.  
    * For each set in the series, computes the average, the variance, its log in base 2.
    */
   // @override
   /*public void testVarianceRate (MonteCarloModelBounded model, int m,
				ArrayList<DensityEstimator> listDE, int numEvalPoints, 
	            double[] integVariance) {
		int n;
		//integVariance = new double[numSets];
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

	   // for(Iterator it=listDE.iterator(); it.hasNext();){
	   // System.out.println(" size"+listDE.size());
	    	for(int i=0; i<listDE.size(); i++) {
	    		System.out.println(" indice "+i);
	  //  for (DensityEstimator de : listDE) {
	 	for (int s = 0; s < numSets; s++) { // For each cardinality n
			n = theSets[s].getNumPoints();
			size[s] = n;
			double[][] data = new double[m][];
			log2n[s] = Num.log2(n);
			// System.out.println(" n = " + n + ", Lg n = " + log2n[s] + "\n"); // ****
			// System.out.println("  " + n + "     " + timer.format());
			RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statReps, data);
			// System.out.println(" data " + Arrays.toString(data[0]));
			//System.out.println(" data1 " + data[0][0]);
			//RQMCExperimentDensity.computeDensityVariance(model, n, m,data, de, numEvalPoints);
			integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model, n, m,
        			data, de, numEvalPoints);
			integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model, n, m,
        			data,listDE.get(i), numEvalPoints);
			
			//RQMCExperimentDensity.computeDensityVarianceListDE (model, n, m, data, listDE, numEvalPoints, integVariance);
			//System.out.println(" integvariance " +integVariance[0]);
			mean[s] = statReps.average();
		    log2Var[s] = Num.log2(integVariance[s]);
		   // System.out.println("  "+log2Var[s]);
		    if (displayExec) {
			   System.out.println("  " + n + "     " + timer.format() + 
			              "   " + PrintfFormat.f(10, 5, mean[s]) + 
					      "   " + PrintfFormat.f(7, 2, log2Var[s]));
		    }
		}	
	 	}   
        cpuTime = timer.format();	   
   }*/
   
   public void testVarianceRate (MonteCarloModelBounded model, int m,
			ArrayList<DensityEstimator> listDE, int numEvalPoints, 
           double[] integVariance, RQMCPointSet [] theSets) {
	   
	   //for (DensityEstimator de : listDE) {
		   for(int i=0; i<listDE.size(); i++) {
		   
		   testVarianceRate (model, m, listDE.get(i), numEvalPoints, integVariance, theSets);
		   
	   }
   
	   
   }
   public void testVarianceRate (MonteCarloModelBounded model, int m,
			DensityEstimator DE, int numEvalPoints, 
           double[] integVariance, RQMCPointSet [] theSets) {
	int n;
	//integVariance = new double[numSets];
	Tally statMean = new Tally();
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

  // for(Iterator it=listDE.iterator(); it.hasNext();){
  // System.out.println(" size"+listDE.size());
   	//for(int i=0; i<listDE.size(); i++) {
   		//System.out.println(" indice "+i);
 //  for (DensityEstimator de : listDE) {
	for (int s = 0; s < numSets; s++) { // For each cardinality n
		n = theSets[s].getNumPoints();
		size[s] = n;
		double[][] data = new double[m][];
		double[][] KS = new double[m][];
		double[][] CVM = new double[m][];
		log2n[s] = Num.log2(n);
		// System.out.println(" n = " + n + ", Lg n = " + log2n[s] + "\n"); // ****
		// System.out.println("  " + n + "     " + timer.format());
		RQMCExperiment.simulReplicatesRQMCSave (model, theSets[s], m, statMean, data);
		// System.out.println(" data " + Arrays.toString(data[0]));
		//System.out.println(" data1 " + data[0][0]);
		//RQMCExperimentDensity.computeDensityVariance(model, n, m,data, de, numEvalPoints);
		/*integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model, n, m,
   			data, de, numEvalPoints);*/
		integVariance[s]=RQMCExperimentDensity.computeDensityVariance (model, m,
   			data, DE, numEvalPoints);
		
		//RQMCExperimentDensity.computeDensityVarianceListDE (model, n, m, data, listDE, numEvalPoints, integVariance);
		//System.out.println(" integvariance " +integVariance[0]);
		mean[s] = statMean.average();
	    log2Var[s] = Num.log2(integVariance[s]);
	   // System.out.println("  "+log2Var[s]);
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
	public String report (boolean details) {
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
		sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");
		sb.append("  Total CPU Time = " + cpuTime + "\n");
		sb.append("-----------------------------------------------------\n");		
		return sb.toString();
	}
	
	/**
	 * Performs an experiment (testVarianceRate) for each point set series in the given list,
	 * and returns a report as a string. 
	 * 
	 * @param model
	 * @param list
	 * @param m
	 */
	/*public String TestRQMCManyPointTypes (MonteCarloModelBounded model, 
			ArrayList<RQMCPointSet[]> list, int m,
			ArrayList<DensityEstimator> listDE, int numEvalPoints, 
            boolean details) {
		StringBuffer sb = new StringBuffer("");
		numReplicates = m;	
		double[] integVariance= new double[numSets];   // Will contain the IV estimates, for each n.
		

		for (int f=0; f<list.size(); f++){
			for (RQMCPointSet[] ptSeries : list) {
				for (int j=0; j<ptSeries.length; j++){
					integVariance= new double[list.size()];
				//integVariance = new double[list.get(f)[j].getPointSet().getNumPoints()];
				init (ptSeries);            
	         	testVarianceRate (model, m, listDE, numEvalPoints, integVariance);
				sb.append (report (details));			
			}}}
		for (RQMCPointSet[] ptSeries : list) {
			init (ptSeries);			
				//integVariance= new double[lng];
				//integVariance =null;
			//integVariance = new double[list.get(f)[j].getPointSet().getNumPoints()];
			            
         	testVarianceRate (model, m, listDE, numEvalPoints, integVariance);
			sb.append (report (details));			
		}
		return sb.toString();
	}*/
	public String TestRQMCManyPointTypes (MonteCarloModelBounded model, 
			ArrayList<RQMCPointSet[]> list, int m,
			ArrayList<DensityEstimator> listDE, int numEvalPoints, 
            boolean details) {
		StringBuffer sb = new StringBuffer("");
		numReplicates = m;	
		double[] integVariance= new double[numSets];   // Will contain the IV estimates, for each n.
		

		/*for (int f=0; f<list.size(); f++){
			for (RQMCPointSet[] ptSeries : list) {
				for (int j=0; j<ptSeries.length; j++){
					integVariance= new double[list.size()];
				//integVariance = new double[list.get(f)[j].getPointSet().getNumPoints()];
				init (ptSeries);            
	         	testVarianceRate (model, m, listDE, numEvalPoints, integVariance);
				sb.append (report (details));			
			}}}*/
		
		for(int i=0; i < listDE.size(); i++) {
		  for (RQMCPointSet[] ptSeries : list) {
			//init (ptSeries);
			
			//theSets=ptSeries;
						
				//integVariance= new double[lng];
				//integVariance =null;
			//integVariance = new double[list.get(f)[j].getPointSet().getNumPoints()];
			            
         	testVarianceRate (model, m, listDE.get(i), numEvalPoints, integVariance, ptSeries);
			sb.append (report (details));			
		  }
	    }
		return sb.toString();
	}
	
	
	/**
	 * Performs an experiment (testVarianceRate) for each point set series in the given list,
	 * and returns a report as a string. 
	 */
	public String toString () {
		return theSets[0].toString();
	}
}