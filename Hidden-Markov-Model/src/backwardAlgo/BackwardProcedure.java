package backwardAlgo;

import hmm_main.HMMEntry;

public class BackwardProcedure {
	
	private static double[][] beta;
	
	// HMM parameters
	double[][] emissionProb;
	double[][] transitionProb;
	double[] initialProb;
	int[] seq;
	
	public BackwardProcedure() {
		beta = new double[HMMEntry.getStates()][HMMEntry.T];
	}
	
	public double[][] getBeta() {
		return beta;
	}
	
	public void computeBeta() {
		
		for(int i = 0; i < HMMEntry.getStates(); i++) {
			beta[i][HMMEntry.T-1] = 0;
		}

		for(int t = HMMEntry.T-2; t >= 0; t--) {
			for(int i = 0; i < HMMEntry.getStates()/2; i++) {
				double sum1 = 0, sum2 = 0;
				for(int j = 0; j < HMMEntry.getStates(); j++) {
					sum1 = sum1 + Math.pow(Math.E, (beta[j][t+1] + transitionProb[i][j] + emissionProb[seq[t+1]][j]));
					sum2 = sum2 + Math.pow(Math.E, (beta[j][t+1] + transitionProb[i + HMMEntry.T][j] + emissionProb[seq[t+1]][j]));
				}
				beta[i][t] = Math.log(sum1);
				beta[i + HMMEntry.T][t] = Math.log(sum2);
			}
		}
	}
	
	public void betaRecurrence(HMMEntry hmmObject, int[] codedSeq) {
		
		emissionProb = hmmObject.getEmissionProb();
		initialProb = hmmObject.getIntialProb();
		transitionProb = hmmObject.getTransitionProb();
		seq = codedSeq;
		
		computeBeta();
	}
}
