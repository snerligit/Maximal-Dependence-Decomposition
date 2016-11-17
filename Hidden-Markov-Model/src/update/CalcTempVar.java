package update;

import backwardAlgo.BackwardProcedure;
import forwardAlgo.ForwardProcedure;
import hmm_main.HMMEntry;

public class CalcTempVar {
	
	double[][] gamma;
	double[][][] eta;
	
	double[][] alpha;
	double[][] beta;
	
	int[] seq;
	
	// HMM parameters
	double[][] emissionProb;
	double[][] transitionProb;
	double[] initialProb;
	
	public CalcTempVar(HMMEntry hmmObject, ForwardProcedure f, BackwardProcedure b, int[] codedSeq) {
		alpha = f.getAlpha();
		beta = b.getBeta();
		
		gamma = new double[HMMEntry.getStates()][HMMEntry.T];
		eta = new double[HMMEntry.getStates()][HMMEntry.getStates()][HMMEntry.T];
		
		emissionProb = hmmObject.getEmissionProb();
		initialProb = hmmObject.getIntialProb();
		transitionProb = hmmObject.getTransitionProb();
		
		seq = codedSeq;
	}
	
	public double[][] getGamma() {
		return gamma;
	}
	
	public double[][][] getEta() {
		return eta;
	}
	
	public double denominator(int t) {
		
		double value = 0;
		for(int j = 0; j < HMMEntry.getStates(); j++) {
			value = value + Math.pow(Math.E, (alpha[j][t] + beta[j][t]));
		}
		return Math.log(value);
	}
	
	public void computeGamma() {
		for(int t = 0; t < HMMEntry.T; t++) {
			double den = denominator(t);
			for(int i = 0; i < HMMEntry.getStates()/2; i++) {
				gamma[i][t] = alpha[i][t] + beta[i][t] - den;
				gamma[i + HMMEntry.T][t] = alpha[i + HMMEntry.T][t] + beta[i + HMMEntry.T][t] - den;
			}
		}
	}
	
	public void computeEta() {
		
		for(int t = 0; t < HMMEntry.T-1; t++) {
			double den = denominator(t);
			for(int i = 0; i < HMMEntry.getStates()/2; i++) {
				for(int j = 0; j < HMMEntry.getStates()/2; j++) {
					eta[i][j][t] = alpha[i][t] + transitionProb[i][j] + beta[j][t+1] + emissionProb[seq[t+1]][j] - den;
					eta[i][j + HMMEntry.T][t] = alpha[i][t] + transitionProb[i][j + HMMEntry.T] + beta[j + HMMEntry.T][t+1] + emissionProb[seq[t+1]][j + HMMEntry.T] - den;
					eta[i + HMMEntry.T][j][t] = alpha[i + HMMEntry.T][t] + transitionProb[i + HMMEntry.T][j] + beta[j][t+1] + emissionProb[seq[t+1]][j] - den;
					eta[i + HMMEntry.T][j + HMMEntry.T][t] = alpha[i + HMMEntry.T][t] + transitionProb[i + HMMEntry.T][j + HMMEntry.T] + beta[j + HMMEntry.T][t+1] + emissionProb[seq[t+1]][j + HMMEntry.T] - den;
				}
			}
		}
	}
}
