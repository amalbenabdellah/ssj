package umontreal.ssj.mcqmctools;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.ArrayList;

import umontreal.ssj.charts.XYLineChart;
import umontreal.ssj.functionfit.LeastSquares;
import umontreal.ssj.gof.GofStat;
import umontreal.ssj.hups.*;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.*;
import umontreal.ssj.stat.density.DEAveragedShiftedHistogram;
import umontreal.ssj.stat.density.DEAveragedShiftedHistogramWeight;
import umontreal.ssj.stat.density.DEHistogram;
import umontreal.ssj.stat.density.DEKernelDensity;
import umontreal.ssj.stat.density.DensityEstimator;
// import umontreal.ssj.stat.list.ListOfTallies;
import umontreal.ssj.util.Chrono;
import umontreal.ssj.util.Num;
import umontreal.ssj.util.PrintfFormat;

/**
 * Provides generic tools to perform RQMC experiments
 * with a simulation model that implements the MCModelDensity interface.
 */

/**
 * @author Pierre L'Ecuyer
 * 
 */
public class RQMCExperimentDensity extends RQMCExperiment {


	/**
	 * Takes data from previous simulation (m replicates, n points each)
	 * and a density estimator de.
	 * Computes and returns an estimate of the integrated variance (IV) 
	 * for this density estimator, obtained by estimating the variance 
	 * at numEvalPoints equidistant points over [a,b] and summing up.
	 */
	public static double computeDensityVariance (MonteCarloModelBounded model, int m,
			double[][] data, DensityEstimator de, int numEvalPoints) {

		double x, y;
		
		double evalPoints[] = new double[numEvalPoints];  // Points at which the density will be evaluated
		double estimDens[] = new double[numEvalPoints];   // Value of the density at those points
		double meanDens[] = new double[numEvalPoints];    // Average value over rep replicates
		double varDens[] = new double[numEvalPoints];     // Variance at each evaluation point
	
		/*Arrays.fill(varDens,0.0);
		Arrays.fill(meanDens,0.0);*/
		for (int i=0; i < numEvalPoints; i++){
			varDens[i] = 0.0;
			meanDens[i] = 0.0;
			
		}
		
		for (int rep = 0; rep < m; rep++) {
			// Estimate the density and evaluate it at eval points
			de.constructDensity(data[rep]);				
			for (int i=0; i < numEvalPoints; i++)
				evalPoints[i] = (de.getB()-de.getA() )/numEvalPoints  ;
			
			de.evalDensity(evalPoints, estimDens, data[rep]);

			
	        // Update the empirical mean and sum of squares of centered observations at each evaluation point.
			for (int j = 0; j < numEvalPoints; j++) { 
				x = estimDens[j];				
				y = x - meanDens[j];
				meanDens[j] += y / (double) (rep+1);		
				varDens[j] += y * (x - meanDens[j]);
				
			}
		}
		double a = model.getMin();
		double b = model.getMax();
		double sumVar = 0.0;
		for (int i = 0; i < numEvalPoints; ++i)
			sumVar += varDens[i];		
		return sumVar * (b - a) / (numEvalPoints * (m - 1));   // Empirical integrated variance.
	}
	
	/*public static double computeDensityBiass (MonteCarloModelDensityKnown model, int m,
			double[][] data, DensityEstimator de, int numEvalPoints) {
		double x,y;
		double meanDens[] = new double[numEvalPoints];
		double biasDens[] = new double[numEvalPoints];
		double evalPoints[] = new double[numEvalPoints];  // Points at which the density will be evaluated
		double estimDens[] = new double[numEvalPoints];   // Value of the density at those points
		
		Arrays.fill(meanDens,0.0);
		//Arrays.fill(biasDens,0.0);

		for (int rep = 0; rep < m; rep++) {
			// Estimate the density and evaluate it at eval points
			de.constructDensity(data[rep]);	
			for (int i=0; i < numEvalPoints; i++)
				evalPoints[i] = (de.getB()-de.getA() )/numEvalPoints  ;
			de.evalDensity(evalPoints, estimDens, data[rep]);

			
	        // Update the empirical mean and sum of squares of centered observations at each evaluation point.
			for (int j = 0; j < numEvalPoints; j++) { 
				x = estimDens[j];				
				y = x - meanDens[j];
				meanDens[j] += y / (double) (rep+1);		
				biasDens[j] = Math.pow(meanDens[j] - model.density(evalPoints[j]),2);
				
			}
			
			}
		double a = model.getMin();
		double b = model.getMax();
		double sumBias = 0.0;
		for (int i = 0; i < numEvalPoints; ++i)
			sumBias += biasDens[i];		
		return sumBias  * (b - a) / (numEvalPoints * (m - 1));   
		
		}*/
	
	
	
	
	/**
	 * Takes data from previous simulation (m replicates, n points each)
	 * and a density estimator de.
	 * Computes and returns an estimate of the mean integrated square error (MISE) 
	 * for this density estimator.
	 */
	
	public static double computeDensityMISE(MonteCarloModelDensityKnown model, int m,
			double[][] data, DensityEstimator de, int numEvalPoints) {

		double x, y;
		
		double evalPoints[] = new double[numEvalPoints];  // Points at which the density will be evaluated
		double estimDens[] = new double[numEvalPoints];   // Value of the density at those points
		double miseDens[] = new double[numEvalPoints];
		//Arrays.fill(miseDens,0.0);
		for (int i=0; i < numEvalPoints; i++){
			miseDens[i] = 0.0;
			
		}
		
		for (int rep = 0; rep < m; rep++) {
			// Estimate the density and evaluate it at eval points
			de.constructDensity(data[rep]);	
			/*if(de==new DEKernelDensity(de.getA(),de.getB())){
			
			for (int i=0; i < numEvalPoints; i++)
				evalPoints[i] = (de.getB()-de.getA() )/numEvalPoints  ;}
			else{*/
				for (int i=0; i < numEvalPoints; i++)
					evalPoints[i] = de.getA() + (i+0.5) *(de.getB()-de.getA() )/numEvalPoints;
			
			de.evalDensity(evalPoints, estimDens, data[rep]);
			
	        // Update the empirical mean and sum of squares of centered observations at each evaluation point.
			for (int j = 0; j < numEvalPoints; j++) { 
				y = estimDens[j] - model.density(evalPoints[j]);   // model.density(x);
				miseDens[j] += y * y;
			}
		}
			
		double a = model.getMin();
		double b = model.getMax();
		
		
		double sumMISE = 0.0;
		for (int i = 0; i < numEvalPoints; ++i)
			sumMISE += miseDens[i];		
		return sumMISE * (b - a) / (numEvalPoints * (m - 1));   
		}
	
	public static double computeDensityBias (MonteCarloModelDensityKnown model, int m,
			double[][] data, DensityEstimator de, int numEvalPoints) {

			
		return   computeDensityMISE (model, m,
				data, de,  numEvalPoints)- computeDensityVariance (model, m,
						data, de,  numEvalPoints);
		}
	

	/*double dx=(model.getMax()-model.getMin())/numEvalPoints;
	evalPoints[0] = min+dx*0.5;
	for (int i=1; i<numEvalPoints; i++){
			evalPoints[i] = evalPoints[i-1]+dx;
			
		}*/
	
	

/*for (int i=0; i < numEvalPoints; i++)
	evalPoints[i] = de.getA() + (i+0.5) *de.geth();
	*/
	/*public static double computeDensityMISEE (MonteCarloModelDensityKnown model, int m,
			double[][] data, DensityEstimator de, int numEvalPoints) {
		
		return   computeDensityVariance (model, m,
				data, de,  numEvalPoints)+ computeDensityBias (model, m,
						data, de,  numEvalPoints);
		}*/
		
		
}
