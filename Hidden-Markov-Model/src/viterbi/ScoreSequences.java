package viterbi;

import hmm_main.HMMEntry;

public class ScoreSequences {
	
	// HMM parameters 
	private double[][] emissionProb;
	private double[][] transitionProb;
	private double[] initialProb;
	
	private int[] Seq;
	private int[] codedSeq;
	private int seqLen;
	
	public ScoreSequences(HMMEntry hmmObject, int[] scoreSeq, int len, int[] coded) {
		emissionProb = hmmObject.getEmissionProb();
		transitionProb = hmmObject.getTransitionProb();
		initialProb = hmmObject.getIntialProb();
		
		Seq = scoreSeq;
		codedSeq = coded;
		seqLen = len;
	}
	
	public double score() {
		double score = initialProb[Seq[0]];
		for(int i = 1; i < seqLen; i++) {
			score += emissionProb[codedSeq[i]][i] + transitionProb[i-1][i];
		}
		return score;
	}

}
