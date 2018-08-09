/*
    GreedyExperimentalDesign
    Software for Experimental Design
    
    Copyright (C) 2015 Professor Adam Kapelner 
    Department of Mathematics, Queens College, City University of New York

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details:
    
    http://www.gnu.org/licenses/gpl-2.0.txt

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

package GreedyExperimentalDesign;

import java.util.ArrayList;

import ExperimentalDesign.*;
import ObjectiveFunctions.*;

/**
 * This class handles initializing many greedy searches for a treatment vector
 * (the design) using a thread pool.
 * 
 * @author Adam Kapelner
 */
public class GreedyExperimentalDesign extends MultipleSearchExperimentalDesigns {
	
	//set by user
	private boolean diagnostics;
	private boolean semigreedy;

	//data inputed from the user's datas
	private Integer max_iters;
	private int[][] starting_indicTs;
	
	//output
	private ArrayList<ArrayList<int[]>> switched_pairs;
	private ArrayList<ArrayList<double[]>> xbardiffjs_by_iterations;
	private ArrayList<ArrayList<Double>> min_obj_val_by_iterations;

	

	//running the Java as standalone is for debug purposes ONLY!!!
	public static void main(String[] args) throws Exception{	

		
		GreedyExperimentalDesign gd = new GreedyExperimentalDesign();
		//set seed here for reproducibility during debugging
		gd.rand_obj.setSeed(1984);

		int n = 300;
		int p = 5;
		gd.setNandP(n, p);
		for (int i = 0; i < n; i++){
//			double[] x_i = {Math.random(), Math.random(), Math.random(), Math.random(), Math.random(), Math.random(), Math.random()};
			double[] x_i = new double[p];
			for (int j = 0; j < p; j++){
				x_i[j] = gd.rand_obj.nextDouble();
			}
			gd.setDataRow(i, x_i);
		}
//		System.out.println("Xstd");
//		for (int i = 0; i < n; i++){
//			System.out.println(Tools.StringJoin(gd.Xstd[i]));
//		}		
		gd.setMaxDesigns(25);
		gd.setObjective(ObjectiveFunction.ABS);
		gd.setDiagnostics();
		gd.setWait();
		gd.beginSearch();
//		System.out.println("progress: " + gd.progress());
	}
	
	public void beginSearch(){
		super.beginSearch();

		//initialize all data
		switched_pairs = new ArrayList<ArrayList<int[]>>(max_designs);
		for (int d = 0; d < max_designs; d++){
			switched_pairs.add(new ArrayList<int[]>());
		}
		xbardiffjs_by_iterations = new ArrayList<ArrayList<double[]>>(max_designs);
		for (int d = 0; d < max_designs; d++){
			xbardiffjs_by_iterations.add(new ArrayList<double[]>());
		}	
		min_obj_val_by_iterations = new ArrayList<ArrayList<Double>>(max_designs);
		for (int d = 0; d < max_designs; d++){
			min_obj_val_by_iterations.add(new ArrayList<Double>());
		}			
//		System.out.println("resulting data initialized");
		
		initializeStartingIndicTs();
		
		for (int d = 0; d < max_designs; d++){
			final int d0 = d;
//			if (d % 100 == 0){
//				System.out.println("worker added to thread pool #" + d);
//			}
	    	search_thread_pool.execute(new Runnable(){
				public void run() {
					new GreedySearch(
							X, 
							Sinv, 
							starting_indicTs[d0], 
							ending_indicTs[d0], 
							switched_pairs.get(d0),
							min_obj_val_by_iterations.get(d0),
							xbardiffjs_by_iterations.get(d0),
							objective_vals, 
							num_iters, 
							objective, 
							d0, 
							semigreedy, 
							diagnostics, 
							max_iters, 
							rand_obj);
				}
			});
		}		
		afterBeginSearch();
		//System.out.println("min_obj_val_by_iterations: " + min_obj_val_by_iterations);
//		System.out.println("num_iters: " + num_iters);
		
	}


	private void initializeStartingIndicTs() {
		starting_indicTs = new int[max_designs][n];
		for (int d = 0; d < max_designs; d++){
			starting_indicTs[d] = Tools.fisherYatesShuffle(Tools.newBalancedBlankDesign(n), rand_obj);
		}
	}	
	
	public int[][][] getSwitchedPairs(int[] indicies){		
		int[][][] pairs_by_iteration_per_search = new int[indicies.length][][];
		for (int i = 0; i < indicies.length; i++){
			ArrayList<int[]> iteration_switched_pairs_raw = switched_pairs.get(indicies[i]);
			int iters = iteration_switched_pairs_raw.size();
			int[][] iteration_switched_pairs = new int[iters][];
			for (int j = 0; j < iters; j++){
				iteration_switched_pairs[j] = iteration_switched_pairs_raw.get(j);
			}
			pairs_by_iteration_per_search[i] = iteration_switched_pairs;
		}
		return pairs_by_iteration_per_search;
	}
	
	public double[][][] getXbarjDiffs(int[] indicies){		
		double[][][] xbarj_diffs_per_search = new double[indicies.length][][];
		for (int i = 0; i < indicies.length; i++){
			ArrayList<double[]> xbarj_diffs_raw = xbardiffjs_by_iterations.get(indicies[i]);
			int iters = xbarj_diffs_raw.size();
			double[][] xbarj_diffs = new double[iters][];
			for (int j = 0; j < iters; j++){
				xbarj_diffs[j] = xbarj_diffs_raw.get(j);
			}
			xbarj_diffs_per_search[i] = xbarj_diffs;
		}
		return xbarj_diffs_per_search;
	}	
	
	public double[][] getObjValByIter(int[] indicies){		
		double[][] min_obj_val_by_iterations_prim = new double[indicies.length][];
		for (int i = 0; i < indicies.length; i++){
			ArrayList<Double> min_obj_val_by_iteration = min_obj_val_by_iterations.get(indicies[i]);
			int iters = min_obj_val_by_iteration.size();
			double[] min_obj_val_by_iteration_prim = new double[iters];
			for (int j = 0; j < iters; j++){
				min_obj_val_by_iteration_prim[j] = min_obj_val_by_iteration.get(j);
			}
			min_obj_val_by_iterations_prim[i] = min_obj_val_by_iteration_prim;
		}
		return min_obj_val_by_iterations_prim;
	}	
	
	public int[][] getStartingIndicTs(int[] indicies){
		int[][] starting_indicTs = new int[indicies.length][n];
		for (int i = 0; i < indicies.length; i++){
			starting_indicTs[i] = this.starting_indicTs[indicies[i]];
		}
		return starting_indicTs;
	}
	
	public void setDiagnostics(){
		diagnostics = true;
	}
	
	public void setSemigreedy(){
		semigreedy = true;
	}
	
	public void setMaxIters(int max_iters){
		this.max_iters = max_iters;
	}
}