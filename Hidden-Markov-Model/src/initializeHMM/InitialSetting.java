package initializeHMM;

import hmm_main.HMMEntry;

public class InitialSetting {
	
	public static int pseudocountN;
	public static int pseudocountD;
	
	public InitialSetting() {
		pseudocountN = 1;
		pseudocountD = 4;
	}
	
	public double[][] formEmissionProbTable(String[] Seq) {
		
		double[][] emissionProb = new double[HMMEntry.getVocabSize()][Seq[0].length()]; 
		
		try {
			for(int k = 0; k < Seq[0].length(); k++) {
				int ctrA = 0, ctrC = 0, ctrG = 0, ctrT = 0;
				for(int j = 0; j < Seq.length; j++) {
					switch(Seq[j].charAt(k)) {
					case 'A': ctrA++; break;
					case 'C': ctrC++; break;
					case 'G': ctrG++; break;
					case 'T': ctrT++; break;
					}
				}
				
				emissionProb[0][k] = Math.log((double)(ctrA + pseudocountN)) - Math.log((Seq.length + pseudocountD));
				emissionProb[1][k] = Math.log((double)(ctrC + pseudocountN)) - Math.log((Seq.length + pseudocountD));
				emissionProb[2][k] = Math.log((double)(ctrG + pseudocountN)) - Math.log((Seq.length + pseudocountD));
				emissionProb[3][k] = Math.log((double)(ctrT + pseudocountN)) - Math.log((Seq.length + pseudocountD));
			}
			
		} catch(Exception e) {
			System.out.println("formEmissionProbTable(): " + e);
		}
		
		return emissionProb;
	}
	
	public double[][] computeEmissionProb(String[] ssSeq, String[] nssSeq) {
		
		double[][] emissionProb = new double[HMMEntry.getVocabSize()][HMMEntry.getStates()];
		double[][] eProbss;
		double[][] eProbnss;
		
		eProbss = formEmissionProbTable(ssSeq);
		eProbnss = formEmissionProbTable(nssSeq);
		
		for(int i = 0; i < HMMEntry.getVocabSize(); i++) {
			for(int j = 0; j < ssSeq[0].length(); j++) {
				emissionProb[i][j] = eProbss[i][j];
				emissionProb[i][j+ssSeq[0].length()] = eProbnss[i][j];
			}
		}
		
		return emissionProb;
	}
	
	public double[][] computeTransitionProb() {
		
		final int DIMENSION = HMMEntry.getStates();
		double[][] transitionProb = new double[DIMENSION][DIMENSION];
		
		//double initialVal = 0 - Math.log((double)DIMENSION);
		double highProb = Math.log(0.8);
		double lowProb = Math.log(0.2);
		
		for(int i = 0; i < DIMENSION; i++) {
			for(int j = 0; j < DIMENSION; j++) {
				if(i == j-1 && i < (DIMENSION/2)) {
					transitionProb[i][j] = highProb;
					transitionProb[i][j+(DIMENSION/2)-1] = lowProb;
				}
				else if(i == j-1 && i >= (DIMENSION/2)) {
					transitionProb[i][j] = highProb;
					transitionProb[i][j-(DIMENSION/2)] = lowProb;
				}
				else {
					transitionProb[i][j] = -99;
				}
			}
		}
		return transitionProb;
	}
	
	public double[] computeInitialProb() {
		
		final int DIMENSION = HMMEntry.getStates();
		
		double[] initialProb = new double[DIMENSION];
		//double initialVal = 0 - Math.log((double)DIMENSION);
		
		initialProb[0] = Math.log(0.2);
		
		for(int i = 1; i < DIMENSION; i++) {
			initialProb[i] = -99;
		}
		
		initialProb[9] = Math.log(0.8);
		return initialProb;
	}
}
