package RerandomizationExperimentalDesign;

import java.util.ArrayList;

import ObjectiveFunctions.*;
import ExperimentalDesign.*;

public class RerandomizationExperimentalDesign extends MultipleSearchExperimentalDesigns {
	
	//set by user
	private Double obj_val_cutoff_to_include;
	
	//running the Java as standalone is for debug purposes ONLY!!!
	public static void main(String[] args) throws Exception{	

		RerandomizationExperimentalDesign rd = new RerandomizationExperimentalDesign();
		//set seed here for reproducibility during debugging
		rd.rand_obj.setSeed(1984);

		int n = 10;
		int p = 20;
		rd.setNandP(n, p);
		for (int i = 0; i < n; i++){
			double[] x_i = new double[p];
			for (int j = 0; j < p; j++){
				x_i[j] = rd.rand_obj.nextDouble();
			}
			rd.setDataRow(i, x_i);
		}
		
		rd.Sinv = new double[p][p];
		for (int i = 0; i < p; i++){
			for (int j = 0; j < p; j++){
				rd.Sinv[i][j] = i == j ? 1 : 0;			
			}
		}
		
		rd.setObjective(ObjectiveFunction.MAHAL);
		rd.setMaxDesigns(1000);
		rd.setObjValCutoffToInclude(10000D);
		rd.setNumCores(3);
		rd.setWait();
		rd.beginSearch();
		
		double[] obj_vals = rd.getObjectiveVals();
		System.out.println("obj_vals: " + Tools.StringJoin(obj_vals));	
		
		for (int i = 0; i < rd.ending_indicTs.length; i++){
			System.out.println("indicT " + (i + 1) + ": " + Tools.StringJoin(rd.ending_indicTs[i]));
		}
	
	}
		
	public void beginSearch(){
		super.beginSearch();
		
		//we gotta calculate the obj function
		ObjectiveFunction obj_fun = null;
		if (objective.equals(ObjectiveFunction.MAHAL)){
			obj_fun = new MahalObjective(Sinv, n);
		}
		else if (objective.equals(ObjectiveFunction.ABS)){
			obj_fun = new AbsSumObjective();	
		}
		
		final ObjectiveFunction fin_obj_fun = obj_fun;

		System.out.println("before pool");

    	search_thread_pool.execute(new Runnable(){
			public void run() {
				while (true){
					int r = progress();
//					System.out.println("progress = " + r);
					//break up here too to avoid one more iteration (ugly, but a tad faster)
					if (r == max_designs){
						break;
					}
					
					int[] indicT = Tools.fisherYatesShuffle(Tools.newBalancedBlankDesign(n), rand_obj);
//					System.out.println("indicT " + Tools.StringJoin(indicT));
					if (obj_val_cutoff_to_include != null){

						int[] i_Ts = Tools.findIndicies(indicT, n / 2, 1);
						int[] i_Cs = Tools.findIndicies(indicT, n / 2, 0);
						ArrayList<double[]> XT = Tools.subsetMatrix(X, i_Ts); 
						ArrayList<double[]> XC = Tools.subsetMatrix(X, i_Cs);
						double[] avg_Ts = Tools.colAvg(XT, p);
						double[] avg_Cs = Tools.colAvg(XC, p);	
						fin_obj_fun.setXTbar(avg_Ts);
						fin_obj_fun.setXCbar(avg_Cs);
						double obj_val = fin_obj_fun.calc(false);
						
						if (obj_val < obj_val_cutoff_to_include){
							if (r == max_designs){
								break;
							}
							//create the new vector and its corresponding objective value
							ending_indicTs[r] = indicT;
							objective_vals[r] = obj_val;
						}
					}
					else {
						//we are just looking for a certain number and then we're done
						if (r == max_designs){
							break;
						}
						ending_indicTs[r] = indicT;
						objective_vals[r] = Double.NaN; //flag it
					}
				}
			}
		});		
		afterBeginSearch();		
	}
	
	public void setObjValCutoffToInclude(double obj_val_cutoff_to_include){
		this.obj_val_cutoff_to_include = obj_val_cutoff_to_include;
	}	
}
