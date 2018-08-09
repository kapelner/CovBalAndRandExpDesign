package ExperimentalDesign;

public abstract class MultipleSearchExperimentalDesigns extends AllExperimentalDesigns {

	protected int max_designs;
	protected int[][] ending_indicTs;	
	protected Double[] objective_vals;	
	protected Integer[] num_iters;
	
	public void beginSearch(){
		super.beginSearch();

		num_iters = new Integer[max_designs];
		ending_indicTs = new int[max_designs][n];
		objective_vals = new Double[max_designs];
	}

	
	public void setMaxDesigns(int max_designs){
		this.max_designs = max_designs;
//		System.out.println("max_designs " + this.max_designs);
	}

	public int progress(){
		int done = 0;
		if (objective_vals != null){
			for (int d = 0; d < max_designs; d++){
//				System.out.println("objective_vals d = " + objective_vals[d]);
				if (objective_vals[d] == null){
					break;
				}
				done++;
			}
		}
		return done;
	}
	
	public int[] getNumIters(){		
		int d_finished = progress();
		int[] num_iters = new int[d_finished];
		for (int i = 0; i < d_finished; i++){
			num_iters[i] = this.num_iters[i];
		}
		return num_iters;
	}
	
	public double[] getObjectiveVals(){		
		int d_finished = progress();
		double[] objective_vals = new double[d_finished];
		for (int i = 0; i < d_finished; i++){
			objective_vals[i] = (this.objective_vals[i] == null) ? 0 : this.objective_vals[i];
		}
		return objective_vals;
	}
	
	public int[][] getEndingIndicTs(int index){ //stupid R
		int[] indicies = {index};
		return getEndingIndicTs(indicies);
	}
	
	public int[][] getEndingIndicTs(int[] indicies){
		int[][] ending_indicTs = new int[indicies.length][n];
		for (int i = 0; i < indicies.length; i++){
			ending_indicTs[i] = this.ending_indicTs[indicies[i]];
		}
		return ending_indicTs;
	}
	
	public int[][] getEndingIndicTs(){
		return ending_indicTs;
	}
}
