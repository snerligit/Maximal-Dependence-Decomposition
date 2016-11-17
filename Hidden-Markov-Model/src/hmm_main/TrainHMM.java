package hmm_main;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import forwardAlgo.ForwardProcedure;
import initializeHMM.InitialSetting;
import update.CalcTempVar;
import update.UpdateParameters;
import viterbi.ViterbiDecoding;
import backwardAlgo.BackwardProcedure;

public class TrainHMM {
	double[][] emissionProb;
	double[][] transitionProb;
	
	int i;
	int merCtr;
	
	public TrainHMM(int i, int merCtr) {
		transitionProb = new double[HMMEntry.getStates()][HMMEntry.getStates()];
		emissionProb = new double[HMMEntry.getVocabSize()][HMMEntry.getStates()];
		this.i = i;
		this.merCtr = merCtr;
	}
	
	public void updateTP(double[][] tpN, double[][] tpD) {
		for(int i = 0; i < HMMEntry.getStates(); i++) {
			for(int j = 0; j < HMMEntry.getStates(); j++) {
				transitionProb[i][j] = Math.log(tpN[i][j]) - Math.log(tpD[i][j]);
			}
		}
	}
	
	public void updateEP(double[][] epN, double[][] epD) {
		for(int i = 0; i < HMMEntry.getVocabSize(); i++) {
			for(int j = 0; j < HMMEntry.getStates(); j++) {
				emissionProb[i][j] = Math.log(epN[i][j]) - Math.log(epD[i][j]);
			}
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
	
	public void startTraining(HMMEntry hmmObject) throws IOException {
		
		InitialSetting randomSet = new InitialSetting();
		hmmObject.setEmissionProb(randomSet.computeEmissionProb(hmmObject.ssSeq, hmmObject.nssSeq));
		hmmObject.setIntialProb(randomSet.computeInitialProb());
		hmmObject.setTransitionProb(randomSet.computeTransitionProb());
	
		System.out.println("\nInitial Prob");
		HMMEntry.printArray(hmmObject.getIntialProb(), HMMEntry.getStates());
		System.out.println("\nEmission Prob");
		HMMEntry.printMatrix(hmmObject.getEmissionProb(), HMMEntry.getVocabSize(), HMMEntry.getStates());
		System.out.println("\nTransition Prob");
		HMMEntry.printMatrix(hmmObject.getTransitionProb(), HMMEntry.getStates(), HMMEntry.getStates());
		
		int iteration;
		
		for(iteration = 0; iteration < HMMEntry.MAX_ITERATIONS; iteration++) {
			
			hmmObject.storeCurrentValues();
		
			for(int steps = 0; steps < hmmObject.ssSeq.length + hmmObject.nssSeq.length; steps++) {
		
				int[] codedSeq;
				if(steps < hmmObject.ssSeq.length) {
					codedSeq = computeCode(hmmObject.ssSeq[steps]);
				}
				else {
					codedSeq = computeCode(hmmObject.nssSeq[steps-hmmObject.ssSeq.length]);
				}
			
				ForwardProcedure f = new ForwardProcedure();
				f.alphaRecurrence(hmmObject, codedSeq);
		
				BackwardProcedure b = new BackwardProcedure();
				b.betaRecurrence(hmmObject, codedSeq);
		
				CalcTempVar c = new CalcTempVar(hmmObject, f, b, codedSeq);
				c.computeGamma();
				c.computeEta();
				
				UpdateParameters u = new UpdateParameters(hmmObject, c, codedSeq, steps);
				u.updateInitialProb();
				u.updateTransitionProb(hmmObject);
				u.updateEmissionProb(hmmObject);
				
			}
			
			updateTP(hmmObject.transitionProbNum, hmmObject.transitionProbDen);
			updateEP(hmmObject.emissionProbNum, hmmObject.emissionProbDen);
			
			if(hmmObject.checkConvergence(transitionProb, emissionProb, iteration)) {
				break;
			}
			
			hmmObject.setEmissionProb(emissionProb);
			hmmObject.setTransitionProb(transitionProb);
			
		}
		
		System.out.println("\nInitial Prob");
		HMMEntry.printArray(hmmObject.getIntialProb(), HMMEntry.getStates());
		System.out.println("\nEmission Prob");
		HMMEntry.printMatrix(hmmObject.getEmissionProb(), HMMEntry.getVocabSize(), HMMEntry.getStates());
		System.out.println("\nTransition Prob");
		HMMEntry.printMatrix(hmmObject.getTransitionProb(), HMMEntry.getStates(), HMMEntry.getStates());
	
		System.out.println("Converging after iterations: " + iteration);
		
		String outFileTrue1 = "/Users/snerli/Work/SJSU/HMM/HBB_Train_Authentic_Test_Cryptic.txt";
		String testFileTrue1 = "/Users/snerli/Work/SJSU/Santrupti/Personal/SJSU/CS298/Datasets/HBB/HBB_5_Positives.txt";
		int trueSite = 1;
		writeToOutput(hmmObject, outFileTrue1, testFileTrue1, trueSite);
		//String outFile = "F:/SJSU/Spring_2015/CS_298/HBB_3_pos.txt";
		//String testFile = "F:/SJSU/Spring_2015/CS_298/DataSets/Testsets/HBB/HBB_3_positives.txt";
		
		// FOR BRCA1
		
		//String outFileTrue = "C:/Santrupti/Personal/SJSU/CS298/Results/HMM/Train_Authentic_Test_HBB/sliding_window_IE_ss.txt";
		//String testFileTrue = "C:/Santrupti/Personal/SJSU/CS298/Datasets/HBB/beta_globin_14_mers.txt";
		//int trueSite = 1;
		
		/* String outFileTrue2 = "C:/Santrupti/Personal/SJSU/CS298/Results/HMM/10-fold-Authentic-Negative-Cryptic-Positive/EI/cbm_ss_" + this.i + ".txt";
		String testFileTrue2 = "C:/Santrupti/Personal/SJSU/CS298/Datasets/cryptic_before_mutation.txt";
		writeToOutput(hmmObject, outFileTrue2, testFileTrue2, trueSite); */
		
		/*String outFileFalse = "/Users/snerli/Desktop/Personal/SJSU/CS298/Results/HMM/EI/" + this.merCtr + "-mer/nss_" + this.i + ".txt";
		String testFileFalse = "/Users/snerli/Desktop/Personal/SJSU/CS298/Datasets/EI/" + this.merCtr + "-mer/EI_negatives_test_" + this.i + ".txt";
		int falseSite = 0;*/
		
		/*String outFileFalse = "C:/Santrupti/Personal/SJSU/CS298/Results/HMM/10-fold-Authentic-Negative-Cryptic-Positive/EI/nss_" + this.i + ".txt";
		String testFileFalse = "C:/Santrupti/Personal/SJSU/CS298/Datasets/EI_true_9.txt";
		int falseSite = 0;
		
		writeToOutput(hmmObject, outFileFalse, testFileFalse, falseSite);*/

	}
	
	public void writeToOutput(HMMEntry hmmObject, String outFile, String testFile, int site) throws IOException {
		
		BufferedWriter brw = new BufferedWriter(new FileWriter(outFile));
		BufferedReader brr = new BufferedReader(new FileReader(testFile));
		
		String line = brr.readLine();
		while(line != null) {
			ViterbiDecoding v = new ViterbiDecoding(hmmObject, line);
			double score = v.decode(hmmObject);
			brw.write(line + "\t" + score + "\t" + site + "\n");
			line = brr.readLine();
			//line = null;
		}
		brw.close();
		brr.close();
	}

}
