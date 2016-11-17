package update;

import hmm_main.HMMEntry;

public class UpdateParameters {
	
	int[] seq;
	
	double[][] gamma;
	double[][][] eta;
	
	int K;
	
	// HMM parameters
	double[] initialProb;
	
	public UpdateParameters(HMMEntry hmmObject, CalcTempVar c, int[] codedSeq, int steps) {
		
		final int DIMENSION = HMMEntry.getStates();
		
		initialProb = new double[DIMENSION];
		
		gamma = c.getGamma();
		eta = c.getEta();
		
		seq = codedSeq;
		
		K = steps;
	}
	
	public void updateInitialProb() {

		for(int i = 0; i < HMMEntry.getStates(); i++) {
			initialProb[i] = gamma[i][1];
		}
	}
	
	public void updateTransitionProb(HMMEntry hmmObject) {
		for(int i = 0; i < HMMEntry.getStates()/2; i++) {
			for(int j = 0; j < HMMEntry.getStates()/2; j++) {
				double etaSum1 = 0, etaSum2 = 0, etaSum3 = 0, etaSum4 = 0;;
				for(int t = 0; t < HMMEntry.T-1; t++) {
					etaSum1 += Math.pow(Math.E, eta[i][j][t]);
					etaSum2 += Math.pow(Math.E, eta[i][j + HMMEntry.T][t]);
					etaSum3 += Math.pow(Math.E, eta[i + HMMEntry.T][j][t]);
					etaSum4 += Math.pow(Math.E, eta[i + HMMEntry.T][j + HMMEntry.T][t]);
				}
				hmmObject.transitionProbNum[i][j] += etaSum1;
				hmmObject.transitionProbNum[i][j + HMMEntry.T] += etaSum2;
				hmmObject.transitionProbNum[i + HMMEntry.T][j] += etaSum3;
				hmmObject.transitionProbNum[i + HMMEntry.T][j + HMMEntry.T] += etaSum4;
				
				double sum1 = sumGamma(i);
				double sum2 = sumGamma(i + HMMEntry.T);
				hmmObject.transitionProbDen[i][j] += sum1;
				hmmObject.transitionProbDen[i][j+ HMMEntry.T] += sum1;
				hmmObject.transitionProbDen[i+ HMMEntry.T][j] += sum2;
				hmmObject.transitionProbDen[i+ HMMEntry.T][j+ HMMEntry.T] += sum2;
			}
		}
	}
	
	public double sumGamma(int i) {
		double sum = 0;
		for(int t = 0; t < HMMEntry.T-1; t++) {
			sum = sum + Math.pow(Math.E, gamma[i][t]);
		}
		return sum;
	}
	
	public void updateEmissionProb(HMMEntry hmmObject) {
		
		for(int k = 0; k < HMMEntry.getVocabSize(); k++) {
			for(int i = 0; i < HMMEntry.getStates()/2; i++) {
				double numerator1 = 0, numerator2 = 0;
				for(int t = 0; t < HMMEntry.T-1; t++) {
					if(seq[t] == k) {
						numerator1 += Math.pow(Math.E, gamma[i][t]);
						numerator2 += Math.pow(Math.E, gamma[i + HMMEntry.T][t]);
					}
				}
				hmmObject.emissionProbNum[k][i] += numerator1;
				hmmObject.emissionProbNum[k][i + HMMEntry.T] += numerator2;
				
				hmmObject.emissionProbDen[k][i] += sumGamma(i);
				hmmObject.emissionProbDen[k][i + HMMEntry.T] += sumGamma(i + HMMEntry.T);
			}
		}
	}
}
