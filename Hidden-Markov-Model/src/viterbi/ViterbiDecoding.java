package viterbi;

import hmm_main.HMMEntry;

public class ViterbiDecoding {
	
	// HMM parameters 
	private double[][] emissionProb;
	private double[][] transitionProb;
	private double[] initialProb;
	
	private String unknownSeq;
	
	int states;
	int[] state;
	
	public ViterbiDecoding(HMMEntry hmmObject, String scoreSeq) {
		emissionProb = hmmObject.getEmissionProb();
		transitionProb = hmmObject.getTransitionProb();
		initialProb = hmmObject.getIntialProb();
		
		states = HMMEntry.getStates();
		unknownSeq = scoreSeq;
		
		state = new int[states];
		for(int i = 0; i < states; i++) {
			state[i] = i;
		}
	}
	
	public int[] computeCode(String seq) {
		
		int[] code = new int[seq.length()];
		
		for(int i = 0; i < seq.length(); i++) {
			switch(seq.charAt(i)) {
			case 'A': code[i] = 0; break;
			case 'C': code[i] = 1; break;
			case 'G': code[i] = 2; break;
			case 'T': code[i] = 3; break;
			}
		}
		return code;
	}
	
	public double decode(HMMEntry hmmObject) {
		
		double[][] T1 = new double[states][unknownSeq.length()];
		int[][] T2 = new int[states][unknownSeq.length()];
		int[] z = new int[unknownSeq.length()];
		int[] x = new int[unknownSeq.length()];
		
		int[] codedSeq = computeCode(unknownSeq);
		
		for(int i = 0; i < states; i++) {
			T1[i][0] = initialProb[i] + emissionProb[codedSeq[0]][i];
			T2[i][0] = 0;
		}
		
		double intermediate;
		for(int i = 1; i < unknownSeq.length(); i++) {
			
			intermediate = 0.0;
			for(int j = 0; j < states; j++)
			{
				T1[j][i] = (T1[0][i-1] + transitionProb[0][j] + emissionProb[codedSeq[i]][j]);
				for(int k = 1; k < states; k++)
				{
					intermediate = (T1[k][i-1] + transitionProb[k][j] + emissionProb[codedSeq[i]][j]);
					if(T1[j][i] < intermediate)
					{
						T1[j][i] = intermediate;
						T2[j][i] = k;
					}
				}
			}
		}
		
		double max = -999.99;
		for(int k = 0; k < states; k++)
		{
			if(max < T1[k][unknownSeq.length()-1])
			{
				max = T1[k][unknownSeq.length()-1];
				z[unknownSeq.length()-1] = k;
			}
		}
		
		x[unknownSeq.length()-1] = state[z[unknownSeq.length()-1]];

		for(int i = unknownSeq.length()-1; i > 0; i--)
		{
			z[i-1] = T2[z[i]][i];
			x[i-1] = state[z[i-1]];
		}
		
		System.out.println("\nDecoded Sequence:");
		for(int i = 0; i < x.length; i++) {
			System.out.print((x[i]+1) + " ");
		}
		System.out.println();
		
		ScoreSequences s = new ScoreSequences(hmmObject, x, unknownSeq.length(), codedSeq);
		double score = s.score();
		return score;
	}
}
