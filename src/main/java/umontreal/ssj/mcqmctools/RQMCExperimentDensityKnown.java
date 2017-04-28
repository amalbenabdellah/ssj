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
import umontreal.ssj.stat.density.DEHistogram;
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
public class RQMCExperimentDensityKnown extends RQMCExperimentDensity {


	/**
	 * Takes data from previous simulation (m replicates, n points each)
	 * and a density estimator de.
	 * Computes and returns an estimate of the integrated variance (IV) 
	 * for this density estimator, obtained by estimating the variance 
	 * at numEvalPoints equidistant points over [a,b] and summing up.
	 */
	
	public static double computeDensityMISE (MonteCarloModelDensityKnown model, int n, int m,
			double[][] data, DensityEstimator de, int numEvalPoints) {

		double x, y;
		// If the density estimator is a histogram, here we may reset numEvalPoints to 
		// the number of bins of the histogram.
		if ( de == new DEHistogram(de.getA(),de.getB()) )
			numEvalPoints = de.getM();
		
		double evalPoints[] = new double[numEvalPoints];  // Points at which the density will be evaluated
		double estimDens[] = new double[numEvalPoints];   // Value of the density at those points
		double meanDens[] = new double[numEvalPoints];    // Average value over rep replicates
		double varDens[] = new double[numEvalPoints];     // Variance at each evaluation point
		double miseDens[] = new double[numEvalPoints];
		
		for (int rep = 0; rep < m; rep++) {
			// Estimate the density and evaluate it at eval points
			de.constructDensity(data[rep]);	
			double min = data[rep][0];
			double max = data[rep][0];
			for (double d : data[rep]) {
				if (d<min) {
					min = d;
				}
				if (d>max) {
					max = d;
				}
			}
			
			/*for (int i=0; i < numEvalPoints; i++)
				evalPoints[i] = de.getA() + (i+1-0.5) * (de.getB()-de.getA())/de.getM();*/
			
			
		
			//double dx=(model.getMax()-model.getMin())/numEvalPoints;
			/*double dx=(max-min)/numEvalPoints;
			evalPoints[0] = min+dx*0.5;
			for (int i=1; i<numEvalPoints; i++){
					evalPoints[i] = evalPoints[i-1]+dx;
					
				}*/
						
			
			for (int i=0; i<numEvalPoints; i++)
			evalPoints[i] = data[rep][i];
			de.evalDensity(evalPoints, estimDens, data[rep]);
			
			/*for (int i=0; i<numEvalPoints; i++)
			
			System.out.println("estimDens"+ estimDens[i]);*/
			
	        // Update the empirical mean and sum of squares of centered observations at each evaluation point.
			for (int j = 0; j < numEvalPoints; j++) { 
				/*x = estimDens[j];				
				y = x - meanDens[j];
				meanDens[j] += y / (double) (rep+1);				
				varDens[j] += y * (x - meanDens[j]);*/	
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
		
}
