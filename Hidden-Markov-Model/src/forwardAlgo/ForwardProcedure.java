package forwardAlgo;

import hmm_main.HMMEntry;

public class ForwardProcedure {
	
	private static double[][] alpha;
	
	// HMM parameters
	double[][] emissionProb;
	double[][] transitionProb;
	double[] initialProb;
	int[] seq;
	
	public ForwardProcedure() {
		alpha = new double[HMMEntry.getStates()][HMMEntry.T];
	}
	
	public double[][] getAlpha() {
		return alpha;
	}
	
	public void computeAlpha() {
		
		for(int i = 0; i < HMMEntry.getStates()/2; i++) {
			alpha[i][0] = initialProb[i] + emissionProb[seq[0]][i];
			alpha[i + HMMEntry.T][0] = initialProb[i + HMMEntry.T] + emissionProb[seq[0]][i + HMMEntry.T];
		}

		for(int t = 0; (t+1) < HMMEntry.T; t++) {
			for(int j = 0; j < HMMEntry.getStates()/2; j++) {
				double evaluate1 = 0.0, evaluate2 = 0.0;
				for(int i = 0; i < HMMEntry.getStates(); i++) {
					evaluate1 = evaluate1 + Math.pow(Math.E, (alpha[i][t] + transitionProb[i][j]));
					evaluate2 = evaluate2 + Math.pow(Math.E, (alpha[i][t] + transitionProb[i][j + HMMEntry.T]));
				}
				alpha[j][t+1] = emissionProb[seq[t+1]][j] + Math.log(evaluate1);
				alpha[j + HMMEntry.T][t+1] = emissionProb[seq[t+1]][j + HMMEntry.T] + Math.log(evaluate2);
			}
		}
	}
	
	public void alphaRecurrence(HMMEntry hmmObject, int[] codedSeq) {
		
		emissionProb = hmmObject.getEmissionProb();
		initialProb = hmmObject.getIntialProb();
		transitionProb = hmmObject.getTransitionProb();
		seq = codedSeq;
		
		computeAlpha();
	}
}
