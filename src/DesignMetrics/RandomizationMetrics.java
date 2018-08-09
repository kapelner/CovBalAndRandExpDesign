package DesignMetrics;

public class RandomizationMetrics {
	//all the vectors
	private int n;
	private int r;
	private int[][] ending_indicTs;
	private double[][] p_hat_same_group;
	private double entropy_metric;
	private double se_metric;
	
	public static void main(String[] args){
		RandomizationMetrics rm = new RandomizationMetrics();

		int[] indicT0 = new int[]{1,0,1,0,1,0,0};
		int[] indicT1 = new int[]{1,1,1,0,1,0,0};
		int[] indicT2 = new int[]{1,0,1,0,1,0,0};
		int[] indicT3 = new int[]{1,0,1,0,1,0,0};
		int[] indicT4 = new int[]{1,0,1,0,1,0,0};

		rm.setNandR(7, 5);
		
		rm.setDesign(0, indicT0);
		rm.setDesign(1, indicT1);
		rm.setDesign(2, indicT2);
		rm.setDesign(3, indicT3);
		rm.setDesign(4, indicT4);
		
//		for (int j = 0; j < rm.r; j++){
//			for (int i = 0; i < rm.n; i++){
//
//				System.out.print(rm.ending_indicTs[j][i] + " ");
//				
//			}
//			System.out.print("\n");
//		}
		
		rm.compute();
	}
	
	public RandomizationMetrics(){}
	
	public void setNandR(int n, int r){
		this.n = n;
		this.r = r;
		ending_indicTs = new int[r][n];
	}
	
	public void setDesigns(int[][] ending_indicTs){
		this.ending_indicTs = ending_indicTs;
	}
	
	public void setDesign(int j0, int[] indicT){
		for (int i = 0; i < n; i++){
			ending_indicTs[j0][i] = indicT[i];
		}			
	}
	
	private void estimatePhats(){
		p_hat_same_group = new double[n][n];
		//for each pair we estimate the probability the randomization
		//produces a different assignment		
		for (int i1 = 0; i1 < n - 1; i1++){
			for (int i2 = i1 + 1; i2 < n; i2++){
				int num_same_group = 0;
				for (int j = 0; j < r; j++){
//					System.out.println("j " + j + " i1 " + i1 + " i2 " + i2 + " val1  " + ending_indicTs[j][i1] + " val2 " + ending_indicTs[j][i2]);
					num_same_group += ((ending_indicTs[j][i1] == ending_indicTs[j][i2]) ? 1 : 0);
				}
				p_hat_same_group[i1][i2] = num_same_group / (double)r;

//				System.out.print("[" + (i1 + 1) + "," + (i2 + 1) + "] ");
//				System.out.print(p_hat_same_group[i1][i2] + " ");
			}
//			System.out.print("\n");
		}
	}
	
	public void compute(){
		estimatePhats();
		
		//this is the probability that a random assignment is the same as another one
		double s_n = (n - 2) / ((double)(2 * n - 2));
		//number of pairs n choose 2 = ...
		int num_pairs = n * (n-1) / 2;
		
		//calculate functions of each p_ij
		double sum_entropies = 0;
		double sum_sqd_dev = 0;
		for (int i1 = 0; i1 < n - 1; i1++){
			for (int i2 = i1 + 1; i2 < n; i2++){
				double p_hat = p_hat_same_group[i1][i2];
				sum_entropies += (probTimesLogProb(p_hat) + probTimesLogProb(1 - p_hat));
				sum_sqd_dev += Math.pow(p_hat - s_n, 2);
			}
		}
		
		//now calculate entropy
		double entropy_norm_factor = s_n * Math.log(s_n) + (1 - s_n) * Math.log(1 - s_n);
		entropy_metric = 1 / ((double) num_pairs) * sum_entropies / entropy_norm_factor;
		
		//now calculate se
		double const_factor = 2 / ((double) n) * Math.sqrt((2 * n - 2) / (double)(n - 2));
		se_metric = const_factor * Math.sqrt(sum_sqd_dev);
	}

	private double probTimesLogProb(double p){
		if (p == 0){ //well known limit convention
			return 0;
		}
		return p * Math.log(p);
	}
	
	public double getRandEntropyMetric() {
		return entropy_metric;
	}

	public double getRandStdErrMetric() {
		return se_metric;
	}

	public double[][] getPhats() {
		return p_hat_same_group;
	}
}
