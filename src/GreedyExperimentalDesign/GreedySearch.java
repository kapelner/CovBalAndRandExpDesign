package GreedyExperimentalDesign;

import java.util.ArrayList;
import java.util.Random;

import ExperimentalDesign.*;
import ObjectiveFunctions.*;

public class GreedySearch {

	private double[][] Xstdscaled;
	private double[][] Xstd;
	private int nT;

	public GreedySearch(double[][] Xstd, 
			double[][] Sinvmat, 
			int[] indicT, 
			int[] ending_indicT, 
			ArrayList<int[]> switched_pairs,
			ArrayList<Double> min_obj_val_by_iteration,
			ArrayList<double[]> xbardiffjs_by_iteration,
			Double[] objective_vals, 
			Integer[] num_iters, 
			String objective, 
			int d0, 
			boolean semigreedy, 
			boolean diagnostics, 
			Integer max_iters, 
			Random r) 
	{
		this.Xstd = Xstd;
		nT = Tools.count(indicT, 1);
		int n = Xstd.length;
		int p = Xstd[0].length;		
		
//		System.out.println("GreedySearch: ready to begin " + d0);
		ObjectiveFunction obj_fun = null;
		if (objective.equals(ObjectiveFunction.MAHAL)){
			obj_fun = new MahalObjective(Sinvmat, n);
		} 
		else if (diagnostics && objective.equals(ObjectiveFunction.ABS)){
			obj_fun = new AbsSumObjectiveWithDiagnostics();	
		}
		else if (objective.equals(ObjectiveFunction.ABS)){
			obj_fun = new AbsSumObjective();	
		}
	
		
		createScaledXstd(); //assume nT = n / 2
//		System.out.println("beginSearch: nT = " + nT + " and nC = " + (n - nT));
		
//		int[] i_Tss = Tools.findIndicies(indicT, nT, 1);
//		System.out.println("i_Ts " + Tools.StringJoin(i_Tss));
		
		Double obj_val = null;		
		
		double min_obj_val = Double.MAX_VALUE;
		
		int iter = 0;
		while (true){
//			System.out.println("iter " + iter);
			iter++;	
//			System.out.println("iter++ " + iter);
			
			int[] indicTmin = null;
			double[] xbardiffjs = null;
//			System.out.println("indicTmin " + indicTmin);
			int[] i_Ts = Tools.findIndicies(indicT, nT, 1);
//			System.out.println("i_Ts " + Tools.StringJoin(i_Ts));
			int[] i_Cs = Tools.findIndicies(indicT, n - nT, 0);
//			System.out.println("i_Cs " + Tools.StringJoin(i_Cs));
			if (semigreedy){ //gotta randomize otherwise inefficient
				i_Ts = Tools.fisherYatesShuffle(i_Ts, r);
				i_Cs = Tools.fisherYatesShuffle(i_Cs, r);
			}
			//build the first avg vectors for speed
			ArrayList<double[]> XT = Tools.subsetMatrix(Xstd, i_Ts); 
			ArrayList<double[]> XC = Tools.subsetMatrix(Xstd, i_Cs); 

			double[] avg_Ts = Tools.colAvg(XT, p);
			double[] avg_Cs = Tools.colAvg(XC, p);
//			System.out.println("INIT XTbar: " + Tools.StringJoin(avg_Ts.getData(), ","));
//			System.out.println("INIT XCbar: " + Tools.StringJoin(avg_Cs.getData(), ","));
			int[] switched_pair = new int[2];
			
//			System.out.println("iter " + iter + " #i_Ts: " + i_Ts.length + " #i_Cs: " + i_Cs.length);
			indices_loop: {
				for (int i_T : i_Ts){
					for (int i_C : i_Cs){
						
						int[] indicT_proposal = indicT.clone();
						//make the single switch
						indicT_proposal[i_T] = 0; //i_T is the new control
						indicT_proposal[i_C] = 1; //i_C is the new treatment
						
						
						/////////////////////OLD INEFFICIENT CODE
//						XT = Tools.subsetMatrix(Xstd, nT, i_Ts, i_T, i_C); 
//						XC = Tools.subsetMatrix(Xstd, nT, i_Cs, i_C, i_T);  
//
//						avg_Ts = Tools.colAvg(XT, p);
//						avg_Cs = Tools.colAvg(XC, p);	
//						
//						obj_fun.setXTbar(avg_Ts);
//						obj_fun.setXCbar(avg_Cs);
						//////////////////////

						updateAvgVec(avg_Ts, i_T, i_C, nT);
						obj_fun.setXTbar(avg_Ts);
						
						updateAvgVec(avg_Cs, i_C, i_T, n - nT);
						obj_fun.setXCbar(avg_Cs);
						
						//calculate our objective function (according to the user's specification)		
						obj_val = obj_fun.calc(false);
						

//						System.out.println("  i_T = " + i_T + " i_C = " + i_C + " obj_val = " + obj_val);
						
						if (obj_val < min_obj_val){
							indicTmin = indicT_proposal;
//							System.out.println("best indicT so far " + Tools.StringJoin(indicTmin));
							min_obj_val = obj_val;
//							System.out.println("switched i_T " + i_T + " and i_C " + i_C);
//							obj_fun.calc(true);
							
//							System.out.println("min_obj_val " + min_obj_val + " for iter " + iter);

							if (diagnostics){
								switched_pair[0] = i_T;
								switched_pair[1] = i_C;
								xbardiffjs = ((AbsSumObjectiveWithDiagnostics)obj_fun).getXbardiffjs().clone();
							}
							
							if (semigreedy){ //semigreedy means as soon as we find improvement, we ditch
								break indices_loop;
							}							
						}
						
						//reset the avg vecs
						updateAvgVec(avg_Ts, i_C, i_T, nT);
						updateAvgVec(avg_Cs, i_T, i_C, n - nT);						
					}	
				}
			}
			//we've finished one iteration by checking every possible switch
			//record this switch only if it is a real switch
			if (diagnostics && indicTmin != null){
				switched_pairs.add(switched_pair);
				xbardiffjs_by_iteration.add(xbardiffjs);
				min_obj_val_by_iteration.add(min_obj_val);
			}
			
			//after searching through every possible switch, we didn't find anything, so break
			if (indicTmin == null){
				break;
			}
			//otherwise - continue and update the binary search vector
			else {
				indicT = indicTmin;
			}

			//we can also be done if we hit our upper limit of iterations
			if (max_iters != null && max_iters == iter){
				break;
			}
		}	
//		System.out.println("after while true");
		
		//search is over; ship back the data now
		for (int i = 0; i < indicT.length; i++){
			ending_indicT[i] = indicT[i];
		}
//		System.out.println("ending_indicT " + Tools.StringJoin(ending_indicT));
		objective_vals[d0] = min_obj_val;
		num_iters[d0] = iter - 1;
//		System.out.println("SEARCH DONE obj_val " + min_obj_val + " iters " + (iter - 1));
	}

	private void createScaledXstd() {
		Xstdscaled = new double[Xstd.length][];
		int p = Xstd[0].length;
		for (int i = 0; i < Xstd.length; i++){
			for (int j = 0; j < p; j++){
				if (Xstdscaled[i] == null){
					Xstdscaled[i] = new double[p];
				}
				Xstdscaled[i][j] = Xstd[i][j] / nT;
			}			
		}
	}

	private void updateAvgVec(double[] avg_vec, int i_remove, int i_add, int nT) {
		double[] obs_to_remove = Xstdscaled[i_remove];
		double[] obs_to_add = Xstdscaled[i_add];
		
		for (int j = 0; j < obs_to_add.length; j++){
			avg_vec[j] = avg_vec[j] - obs_to_remove[j] + obs_to_add[j];
		}		
	}
}
