package hmm_main;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class HMMEntry {

	public String[] ssSeq;
	public String[] nssSeq;
	public String method;
	private static int noSsSeq;
	private static int noNssSeq;
	private static int states;
	private static int vocabSize;
	public static int[] codeSeq = {0, 1, 2, 3};
	public static int T;
	public static final int MAX_ITERATIONS = 20;
	
	
	// HMM Parameters
	private double[][] emissionProb;
	private double[][] transitionProb;
	private double[] initialProb;
	
	// Old HMM Parameters
	private double[][] emissionProbOld;
	private double[][] transitionProbOld;
	
	//Keep track of convergence
	double transitionDiff;
	double emissionDiff;
	
	// For updating multiple sequence observations
	public double[][] emissionProbNum;
	public double[][] emissionProbDen;
	public double[][] transitionProbNum;
	public double[][] transitionProbDen;
	
	public HMMEntry() {
		ssSeq = null;
		method = "HMM";
		noSsSeq = 0;
		noNssSeq = 0;
		vocabSize = 4;
		transitionDiff = 0.0;
		emissionDiff = 0.0;
	}
	
	public void setSsNssSeqNo(String ssFile, String nssFile) throws FileNotFoundException, IOException {
		// open the training file to count the length of each sequence and total number of sequences
	
		BufferedReader brSs = new BufferedReader(new FileReader(ssFile));
		BufferedReader brNss = new BufferedReader(new FileReader(nssFile));
		
		try{
			String line = brSs.readLine();
			int ctr = 0;
			
			T = line.length();
			
			while(line != null) {
				if(line.charAt(0) != '>') {
					ctr++;
				}
				line = brSs.readLine();
			}
			
			noSsSeq = ctr;
			
			ctr = 0;
			line = brNss.readLine();
			
			while(line != null) {
				if(line.charAt(0) != '>') {
					ctr++;
				}
				line = brNss.readLine();
			}
			
			noNssSeq = ctr;
			
		} catch(Exception e){
			System.out.println("setSsSeqNo(): " + e);
		} finally {
			brSs.close();
			brNss.close();
		}
	}
	
	public int getTotalSeq() {
		return (noSsSeq + + noNssSeq);
	}
	
public void ReadSsNssFile(String ssFile, String nssFile) throws FileNotFoundException, IOException  {
		
		// Read the splice site set into ssSeq and non-splice site set into nssSeq
		
		BufferedReader brSs = new BufferedReader(new FileReader(ssFile));
		BufferedReader brNss = new BufferedReader(new FileReader(nssFile));
		
		try{
			String line;
			int ctr = 0;
			
			line = brSs.readLine(); 
			
			setStates(line.length() * 2);
			System.out.println("No. of states: " + states);
			
			ssSeq = new String[noSsSeq];
			
			while(line != null) {
				ssSeq[ctr] = new String(line.toUpperCase());
				ctr++;
				line = brSs.readLine();
			}
			
			line = brNss.readLine(); 
			ctr = 0;
			nssSeq = new String[noNssSeq];
			
			while(line != null) {
				nssSeq[ctr] = new String(line.toUpperCase());
				ctr++;
				line = brNss.readLine();
			}
	
		} catch(Exception e){
			System.out.println("ReadSsNssFile(): " + e);
		} finally {
			brSs.close();
			brNss.close();
		}
	}
	
	public void parseCliArgs(int i, int merCtr) {
		try {
			/*if(args.length < 2 || args[0].equalsIgnoreCase("-h")) {
				System.out.println("Usage: java <class_name> <splice_site_file_path> <non_splice_site_file_path>");
				System.exit(0);
			}
			String ssFilePath = args[0];
			String nssFilePath = args[1];*/
			
			//String ssFilePath = "C:/Santrupti/Personal/SJSU/CS298/Datasets/EI/9-mer/EI_positives_train_" + i + ".txt";
			//String nssFilePath = "C:/Santrupti/Personal/SJSU/CS298/Datasets/EI/9-mer/EI_negatives_train_" + i + ".txt";
			
			String ssFilePath = "/Users/snerli/Work/SJSU/Santrupti/Personal/SJSU/CS298/Datasets/EI_true_9.txt";
			String nssFilePath = "/Users/snerli/Work/SJSU/Santrupti/Personal/SJSU/CS298/Datasets/EI_false_9.txt";
			
			//String ssFilePath = "F:/SJSU/Spring_2015/CS_298/DataSets/Testsets/HBB/HBB_3_positives.txt";
			//String nssFilePath = "F:/SJSU/Spring_2015/CS_298/DataSets/Testsets/HBB/HBB_3_negatives.txt";

			System.out.println("Splice Site File - " + ssFilePath);
			System.out.println("Non-splice Site File - " + nssFilePath);
		
			this.setSsNssSeqNo(ssFilePath, nssFilePath);
			this.ReadSsNssFile(ssFilePath, nssFilePath);
			
			emissionProbNum = new double[getVocabSize()][getStates()];
			emissionProbDen = new double[getVocabSize()][getStates()];
			transitionProbNum = new double[getStates()][getStates()];
			transitionProbDen = new double[getStates()][getStates()];
			
			emissionProbOld = new double[HMMEntry.getVocabSize()][HMMEntry.getStates()];
			transitionProbOld = new double[HMMEntry.getStates()][HMMEntry.getStates()];
			
		} catch (Exception e) {
			System.out.println("parseArgs(): " + e);
		}
	}
	
	public void printSeq() {
		for(int i = 0; i < ssSeq.length; i++) {
			System.out.println(ssSeq[i]);
		}
		
		System.out.println("Splice sites end.");
	}
	
	public void setEmissionProb(double[][] matrix) {
		emissionProb = matrix;
	}
	
	public double[][] getEmissionProb() {
		return emissionProb;
	}
	
	public void setTransitionProb(double[][] matrix) {
		transitionProb = matrix;
	}
	
	public double[][] getTransitionProb() {
		return transitionProb;
	}
	
	public void setIntialProb(double[] matrix) {
		initialProb = matrix;
	}
	
	public double[] getIntialProb() {
		return initialProb;
	}
	
	public void setStates(int s) {
		states = s;
	}
	
	public static int getStates() {
		return states;
	}
	
	public static int getVocabSize() {
		return vocabSize;
	}
	
	public static int getNoSeq() {
		return noSsSeq;
	}
	
	public void storeCurrentValues() {
		
		for(int i = 0; i < getStates(); i++) {
			for(int j = 0; j < getStates(); j++) {
				transitionProbOld[i][j] = transitionProb[i][j]; 
			}
		}
		
		for(int i = 0; i < HMMEntry.getVocabSize(); i++) {
			for(int j = 0; j < getStates(); j++) {
				emissionProbOld[i][j] = emissionProb[i][j]; 
			}
		}
	}
	
	public boolean checkConvergence(double[][] tP, double[][] eP, int iteration) {
		double transitionDiff = 0.0;
		double emissionDiff = 0.0;
		for(int i = 0; i < getStates(); i++) {
			for(int j = 0; j < getStates(); j++) {
				transitionDiff += Math.abs(tP[i][j]-transitionProbOld[i][j]); 
			}
		}
		
		for(int i = 0; i < HMMEntry.getVocabSize(); i++) {
			for(int j = 0; j < getStates(); j++) {
				emissionDiff += Math.abs(eP[i][j]-emissionProbOld[i][j]); 
			}
		}
		
		if(iteration != 0 && (transitionDiff > this.transitionDiff || emissionDiff > this.emissionDiff)) {
			return true;
		}
		
		this.transitionDiff = transitionDiff;
		this.emissionDiff = emissionDiff;
		
		return false;
	}
	
	public static void printMatrix(double[][] matrix, int row, int col) {
		for(int i = 1; i <= col; i++) {
			System.out.print("\t" + i + "\t");
		}
		System.out.println();
		
		for(int i = 0; i < row; i++) {
			System.out.print(i+1 + "\t");
			for(int j = 0; j < col; j++) {
				System.out.printf("%f\t", matrix[i][j]);
			}
			System.out.println();
		}
	}
	
	public static void printArray(double[] matrix, int row) {		
		for(int i = 0; i < row; i++) {
				System.out.printf("%f\t", matrix[i]);
		}
	}
	
	public void setHMMParameters(double[][] eP, double[][] tP, double[] iP) {
		this.emissionProb = eP;
		this.transitionProb = tP;
		this.initialProb = iP;
	}
	
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub

		//for (int merCtr = 7; merCtr <= 17; merCtr+=2) {
			for(int i = 1; i <= 1; i++) {
				HMMEntry hmmObject = new HMMEntry();
				//hmmObject.parseCliArgs(args);
				//hmmObject.parseCliArgs(0);
				//hmmObject.parseCliArgs(i, merCtr);
				hmmObject.parseCliArgs(i, 0);
		
				//TrainHMM train = new TrainHMM(0);
				//TrainHMM train = new TrainHMM(i, merCtr);
				TrainHMM train = new TrainHMM(i, 0);
				train.startTraining(hmmObject);
			}
		//}
	}

}
