package umontreal.ssj.mcqmctools;
import umontreal.ssj.gof.GofStat;
import umontreal.ssj.hups.*;
import umontreal.ssj.stat.*;
import umontreal.ssj.stat.list.lincv.ListOfTalliesWithCV;

/**
 * Provides basic generic tools to perform RQMC experiments
 * with a simulation model that implements the MonteCarloModelDouble interface.
 */

/**
 * @author Pierre L'Ecuyer
 * 
 */
public class RQMCExperiment extends MonteCarloExperiment {


	/**
	 * Simulate <tt>m/tt> replications with prqmc and return the m RQMC average observations 
	 * in statReps.
	 */
	public static void simulReplicatesRQMC (MonteCarloModelDouble model, RQMCPointSet prqmc, 
	        int m, Tally statReps) {
		statReps.init();
		int n = prqmc.getNumPoints();
		// Internal collector for stats on the n outputs X, for each replication.
		Tally statValue = new Tally();
		PointSetIterator stream = prqmc.iterator();
		for (int rep = 0; rep < m; rep++) {
			prqmc.randomize();
			stream.resetStartStream();
			simulateRuns(model, n, stream, statValue);
			statReps.add(statValue.average());   // For the estimator of the mean.
		}
	}

	/**
	 * Simulate <tt>m/tt> replications and return the m RQMC average observations in statReps.
	 */
	public static void simulReplicatesRQMC (MonteCarloModelDouble model, PointSet p, PointSetRandomization rand,
	        int m, Tally statReps) {
		statReps.init();
		int n = p.getNumPoints();
		// Internal collector for stats on the n outputs X, for each replication.
		Tally statValue = new Tally();
		PointSetIterator stream = p.iterator();
		for (int rep = 0; rep < m; rep++) {
			rand.randomize(p);
			stream.resetStartStream();
			simulateRuns(model, n, stream, statValue);
			statReps.add(statValue.average());   // For the estimator of the mean.
		}
	}
	
	
		public static void mapRunsTo01 (MonteCarloModelDensityU01 model, 
				TallyStore statValue, TallyStore statValue01) {
			statValue01.init();
			for (int i = 0; i < statValue.numberObs(); i++) 
				statValue01.add(model.cdf(statValue.getArray()[i]));
		}
		
		public static double computeKS(double[] x) {
			double[] ks = GofStat.kolmogorovSmirnov(x);
			return ks[2];
		}

		public static double computeSquareL2StarDim1(double[] x) {
			int n = x.length;
			double discL2Star = 0;
			double[] w = new double[n];
			for (int j = 0; j < n; ++j) {
				w[j] = (j + 0.5) / n;
			}
			// Arrays.sort(x);
			for (int j = 0; j < n; ++j) {
				discL2Star += (w[j] - x[j]) * (w[j] - x[j]);
			}
			discL2Star /= n;
			discL2Star += 1.0 / (12.0 * n * n);
			return discL2Star;
		}


	/**
	 * Simulate <tt>m/tt> replications with prqmc and returns all observations of each
	 * replication in a two-dimensional array: an array of <tt>m/tt> arrays of size <tt>n/tt>.
	 * Each array of size <tt>n/tt> is sorted by increasing order.
	 * This is useful for density and cdf estimation.
	 */
		
		public static void simulReplicatesRQMCSaveU01 (MonteCarloModelDensityU01 model, RQMCPointSet prqmc, 
		        int m, Tally statReps, TallyStore statKS, TallyStore statCVM, double[][] data) {
			int n = prqmc.getNumPoints();
			// Internal collector for stats on the n outputs X, for each replication.
			
			TallyStore statSave = new TallyStore(n);
			TallyStore statSave01 = new TallyStore(n);
			PointSetIterator stream = prqmc.iterator();
			for (int rep = 0; rep < m; rep++) {			
				
				/*if(prqmc.getPointSet() instanceof DigitalNet || prqmc.getPointSet() instanceof CachedPointSet)
				prqmc.randomize();
				else{
				  PointSetRandomization rand = prqmc.getRandomization();
				  prqmc.getPointSet().randomize(rand);
				}*/
				prqmc.randomize();
				
				stream.resetStartStream();			
				simulateRuns(model, n, stream, statSave);
				statSave.quickSort();
				statReps.add(statSave.average());   // For the estimator of the mean.
				data[rep] = statSave.getArray();   // Instead of copy, just exchange pointers!
				
				mapRunsTo01(model,statSave, statSave01);
				double ks=computeKS(statSave01.getArray());
				statKS.add(ks*ks);
				statCVM.add(computeSquareL2StarDim1(statSave01.getArray()));
				statSave = new TallyStore(n);
			    statSave01 = new TallyStore(n);
			    statKS = new TallyStore();
				statCVM = new TallyStore();
			}
		}
		
		/** Performs n runs using stream and collects statistics in statValue.
			Note that the stream can be an RQMC stream.
			*/
		
	public static void simulReplicatesRQMCSave (MonteCarloModelBounded model, RQMCPointSet prqmc, 
	        int m, Tally statMean, double[][] data ) {
		int n = prqmc.getNumPoints();
		// Internal collector for stats on the n outputs X, for each replication.
	
     	TallyStore statSave = new TallyStore(n);
		PointSetIterator stream = prqmc.iterator();
		for (int rep = 0; rep < m; rep++) {			
			
			/*if(prqmc.getPointSet() instanceof DigitalNet || prqmc.getPointSet() instanceof CachedPointSet)
			prqmc.randomize();
			else{
			  PointSetRandomization rand = prqmc.getRandomization();
			  prqmc.getPointSet().randomize(rand);
			}
			*/
			prqmc.randomize();
			stream.resetStartStream();			
			simulateRuns(model, n, stream, statSave);
			statSave.quickSort();
			statMean.add(statSave.average());   // For the estimator of the mean.			
			data[rep] = statSave.getArray();   // Instead of copy, just exchange pointers!			
			statSave = new TallyStore(n);
		}
	}

	/**
	 * Makes <SPAN CLASS="MATH"><I>m</I></SPAN> independent randomizations of the point set
	 * <SPAN CLASS="MATH"><I>p</I></SPAN> using stream <TT>noise</TT>. For each of them, performs
	 * one simulation run for each point of <SPAN CLASS="MATH"><I>p</I></SPAN>, and adds the
	 * averages of the performance (value) and of the CVs over these points to 
	 * <TT>statWithCV</TT>.  This implements the replication method for RQMC with CVs,
	 * as in Section 5 of Hickernell, Lemieux and Owen (2005).
	 * CV estimates, variances, and covariances can be recovered from  <TT>statWithCV</TT>.
	 */
	public static void simulReplicatesRQMCCV (MonteCarloModelCV model, RQMCPointSet prqmc, 
		    int m, ListOfTalliesWithCV<Tally> statWithCV) {
		// Makes m independent randomizations of the point set p using stream
		// noise. For each of them, performs one simulation run for each point
		// of p, and adds the averages over these points to the statWithCV.

		statWithCV.init();
		int n = prqmc.getNumPoints();
		int numCV = model.getNumberCV();
		PointSetIterator stream = prqmc.iterator();
		for (int rep = 0; rep < m; rep++) {
			prqmc.randomize();
			stream.resetStartStream();
			double sumValues = 0.0;
			double[] sumValuesCV = new double[numCV]; 
			double[] curValuesCV = new double[numCV];  // Current value of the CV
			for (int i = 0; i < n; i++) {
				model.simulate(stream);
				sumValues += model.getPerformance();
				curValuesCV = model.getValuesCV();
				for (int k=0; k < numCV; k++) sumValuesCV[k] += curValuesCV[k];
				stream.resetNextSubstream();
			}
			for (int k=0; k < numCV; k++) 
				sumValuesCV[k] /= (double)n;
			statWithCV.add (sumValues/(double)n, sumValuesCV);
		}
	}
}
