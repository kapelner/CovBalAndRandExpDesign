/*
    OptimalExperimentalDesign
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

package OptimalExperimentalDesign;

import ObjectiveFunctions.*;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;

import ExperimentalDesign.AllExperimentalDesigns;
import ExperimentalDesign.Tools;

/**
 * This class handles initializing many optimal searches for a treatment vector
 * (the design) using a thread pool.
 * 
 * @author Adam Kapelner
 */
public class OptimalExperimentalDesign extends AllExperimentalDesigns {

	private static final double NOT_REACHED_YET = -999999;

	private static final int BATCH_SIZE = 100000;
	
	//temp stuff
	private static HashMap<Integer, ArrayList<BitSet>> all_indicTs;
	static {
		all_indicTs = new HashMap<Integer, ArrayList<BitSet>>();
	}
	private int max_designs;
	private int n_over_two;
	
	//output
	private double[] objective_vals;
	private int d_opt;

	//running the Java as standalone is for debug purposes ONLY!!!
	public static void main(String[] args) throws Exception{	

		for (int g = 1; g <= 1; g++){
			OptimalExperimentalDesign od = new OptimalExperimentalDesign();
			//set seed here for reproducibility during debugging
			od.rand_obj.setSeed(1984);
	
			int n = 6;
			int p = g;
			od.setNandP(n, p);
			for (int i = 0; i < n; i++){
	//			double[] x_i = {Math.random(), Math.random(), Math.random(), Math.random(), Math.random(), Math.random(), Math.random()};
				double[] x_i = new double[p];
				for (int j = 0; j < p; j++){
					x_i[j] = od.rand_obj.nextDouble();
				}
				od.setDataRow(i, x_i);
			}
//			System.out.println("Xstd");
//			for (int i = 0; i < n; i++){
//				System.out.println(Tools.StringJoin(od.X[i]));
//			}
			od.setObjective(ObjectiveFunction.ABS);
			od.setNumCores(3);
			od.setWait();
			od.beginSearch();
		}
//		System.out.println("progress: " + od.progress());
	}
	
	public void beginSearch(){
		super.beginSearch();
		
		initializeStartingIndicTs();
		
		objective_vals = new double[max_designs];
		for (int d = 0; d < max_designs; d++){
			objective_vals[d] = NOT_REACHED_YET;
		}
		
		for (int d = 0; d < max_designs; d += BATCH_SIZE){
			final int d0 = d;

	    	search_thread_pool.execute(new Runnable(){
				public void run() {
//					System.out.println("RUN" + all_indicTs.get(n).size());
					for (int d00 = d0; d00 < d0 + BATCH_SIZE; d00++){
						if (d00 >= max_designs){
							break;
						}
//						System.out.println("d00 " + d00 + " stop " + stop);
						if (d00 % 1000000 == 0){
							System.out.println("million");
						}
						ObjectiveFunction obj_fun = null;
						if (objective.equals(ObjectiveFunction.MAHAL)){
							obj_fun = new MahalObjective(Sinv, n);
						}
						else if (objective.equals(ObjectiveFunction.ABS)){
							obj_fun = new AbsSumObjective();	
						}
						
						//get the vector for this run
						BitSet indicTbit = all_indicTs.get(n).get(d00);
//						System.out.println((d0 + 1) + " bitvector: " + Tools.StringJoin(indicTbit, ""));
						int[] indicT = Tools.convert_bitvector_to_intvector(indicTbit, n);
//						System.out.println((d0 + 1) + "intvector: " + Tools.StringJoin(indicT, ""));
						//get the indicies
						int[] i_Ts = Tools.findIndicies(indicT, n_over_two, 1);
//						System.out.println("i_Ts " + Tools.StringJoin(i_Ts));
						int[] i_Cs = Tools.findIndicies(indicT, n - n_over_two, 0);
						//get the rows for each group
						ArrayList<double[]> XT = Tools.subsetMatrix(X, i_Ts); 
						ArrayList<double[]> XC = Tools.subsetMatrix(X, i_Cs); 
						//compute the averages
						double[] avg_Ts = Tools.colAvg(XT, p);
//						System.out.println("avg_Ts " + avg_Ts[0]);
						double[] avg_Cs = Tools.colAvg(XC, p);
//						System.out.println("avg_Cs " + Tools.StringJoin(avg_Cs));
						obj_fun.setXTbar(avg_Ts);
						obj_fun.setXCbar(avg_Cs);
						
						//calculate our objective function (according to the user's specification)
						if (d00 < max_designs){
							objective_vals[d00] = obj_fun.calc(false);
						}
//						System.out.println("objval: " + objective_vals[d00]);
						//break out if user desires
						if (search_stopped){
							break;
						}						
					}
				}
			});
		}
		afterBeginSearch();
		
//		for (int i = 0; i < max_designs; i++){
//			System.out.println((i + 1) + ": " + objective_vals[i]);
//		}
		//now return the min
		d_opt = Tools.min_index(objective_vals);
		System.out.println("size of space: " + max_designs);
		System.out.println("d_opt: " + (d_opt + 1));
		System.out.println("obj_opt: " + objective_vals[d_opt]);
		System.out.println("time elapsed in sec: " + timeElapsedInSeconds());
	}
	
	private void initializeStartingIndicTs() {
		if (all_indicTs.containsKey(n)){
			max_designs = all_indicTs.get(n).size();
			return;
		}		
//		System.out.println("begin initializeStartingIndicTs");
		max_designs = (int)n_choose_k(n, n / 2);
//		System.out.println("max_designs " + max_designs);

		all_indicTs.put(n, new ArrayList<BitSet>(max_designs));

		recursivelyFindAllBinaryVecs(new BitSet(), 0, 0, 0);
//		for (int i = 0; i < max_designs; i++){
//			System.out.println((i + 1) + ": " + Tools.StringJoin(all_indicTs.get(n).get(i), ""));
//		}
//		System.out.println("end initializeStartingIndicTs");
	}

	private void recursivelyFindAllBinaryVecs(BitSet bitSet, int pos, int on, int off) {
		
		//if we've made it to the end, we're done
		if (pos == n){
			all_indicTs.get(n).add(bitSet);
//			System.out.println("bitSet: " + Tools.StringJoin(Tools.convert_bitvector_to_intvector(bitSet, n)));
			if (all_indicTs.size() % 1000000 == 0){
				System.out.println("million");
			}
//			System.out.println(Tools.StringJoin(bitSet, ""));
			return;
		}
		
		//now set the next position on and recurse
		BitSet next_vec_on = (BitSet)bitSet.clone();
		next_vec_on.set(pos, true);
		if (on + 1 <= n_over_two){
			recursivelyFindAllBinaryVecs(next_vec_on, pos + 1, on + 1, off);
		}		
		
		//now set the next position off and recurse
		BitSet next_vec_off =  (BitSet)bitSet.clone();
		next_vec_off.set(pos, false);
		if (off + 1 <= n_over_two){
			recursivelyFindAllBinaryVecs(next_vec_off, pos + 1, on, off + 1);
		}
	}

	//these two methods save me having to import the Apache math library
	private long n_choose_k(int n, int k) {		
		 long n_choose_k = (long)Math.floor(Math.exp(ln_factorial(n) - ln_factorial(k) - ln_factorial(n - k)));
//		 System.out.println(n + " choose " + k + " = " + n_choose_k);
		 return n_choose_k;
	}

	private double ln_factorial(int n) {
		double sum = 0;
		for (int i = 1; i <= n; i++){
			sum += Math.log(i);
		}
//		System.out.println("ln " + n + "! = " + sum);
		return sum;
	}

	private int num_vectors_checked(){
//		System.out.println("max_designs " + max_designs);
		int done = 0;
		if (objective_vals != null){
			for (int d = 0; d < max_designs; d++){
//				System.out.println("progress loop d = " + d);
				if (objective_vals[d] == NOT_REACHED_YET){
					break;
				}
				done++;
			}
		}
		return done;
	}
	
	public double progress(){
		return num_vectors_checked() / (double)objective_vals.length;
	}	
	
	public double[] getObjectiveVals(){		
//		int d_finished = num_vectors_checked();
//		double[] objective_vals = new double[d_finished];
//		for (int i = 0; i < d_finished; i++){
//			objective_vals[i] = this.objective_vals[i];
//		}
		return this.objective_vals;
	}
	public double getOptObjectiveVal(){		
		return objective_vals[d_opt];
	}	
	
	public double[] getAllObjectiveVals(){		
		return objective_vals;
	}
	
	public int[] getOptIndicT() {
		return Tools.convert_bitvector_to_intvector(all_indicTs.get(n).get(d_opt), n);
	}
	
	public void setNandP(int n, int p) throws Exception {
		super.setNandP(n, p);
		n_over_two = n / 2;
	}

	

}