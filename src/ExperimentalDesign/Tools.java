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

package ExperimentalDesign;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collection;
import java.util.Random;

/**
 * A class that contains many generally useful convenience methods.
 * 
 * @author Adam Kapelner
 */
public class Tools {	

	public static double[] colAvg(ArrayList<double[]> X, int p) {
		int n = X.size();	
		double[] tally = new double[p];
		for (int i = 0; i < n; i++){
			for (int j = 0; j < p; j++){
				tally[j] += X.get(i)[j];
			}			
		}
		for (int j = 0; j < p; j++){
			tally[j] /= n;
		}
		return tally;
	}

	public static ArrayList<double[]> subsetMatrix(double[][] X, int[] indices) {
//		System.out.println("subsetMatrix indices: " + StringJoin(indices));
		ArrayList<double[]> Xsub = new ArrayList<double[]>(indices.length);
		for (int i : indices){
			Xsub.add(X[i]);
//			System.out.println(StringJoin(X[i]));
		}
		return Xsub;
	}
	
	public static ArrayList<double[]> subsetMatrix(double[][] Xstd, int nT, int[] indices, int i_remove, int i_add) {
		ArrayList<double[]> XstdT = new ArrayList<double[]>(nT);
		for (int i : indices){
			if (i != i_remove){
				XstdT.add(Xstd[i]);
			}			
		}
		XstdT.add(Xstd[i_add]);		
		return XstdT;
	}

	public static int count(int[] indicT, int val) {
		int tally = 0;
		for (int i = 0; i < indicT.length; i++){
			if (indicT[i] == val){
				tally++;
			}			
		}
		return tally;
	}

	public static int[] findIndicies(int[] vec, int n_val, int val) {
//		System.out.println("findIndicies veclen = " + vec.length + " n_val = " + n_val + " val = " + val);
//		System.out.println(StringJoin(vec));
		int[] indicies = new int[n_val];
		int index = 0;
		for (int i = 0; i < vec.length; i++){			
			if (vec[i] == val){				
				indicies[index] = i;
				index++;
//				System.out.println("    tot = " + index);
			}
//			System.out.println("i = " + i + " done");
		}
//		System.out.println("  indicies " + Tools.StringJoin(indicies));
		return indicies;
	}	

	
	//from http://algs4.cs.princeton.edu/11model/Knuth.java.html
	public static int[] fisherYatesShuffle(int[] arr, Random rand){
	    int n = arr.length;
        for (int i = 0; i < n; i++) {
            // choose index uniformly in [i, N-1]
            int r = i + (int) (rand.nextDouble() * (n - i));
            int swap = arr[r];
            arr[r] = arr[i];
            arr[i] = swap;
        }
	    return arr;
	}
	
	public static int[] newBalancedBlankDesign(int n){
		int[] design = new int[n];
		for (int i = 0; i < n; i++){
			design[i] = i < n / 2 ? 1 : 0;
		}
		return design;
	}
	
	/**
	 * Joins a collection of strings into one string
	 * 
	 * @param all		the collection of substrings
	 * @param joinby	the token that joins the substrings
	 * @return			the final product: str1 + joinby + str2 + . . . + strN
	 */	
	@SuppressWarnings("rawtypes")
	public static String StringJoin(ArrayList all, String joinby){
		if (all == null){
			return " NULL ARRAY ";
		}		
		return StringJoin(all.toArray(), joinby);
	}	
	
	/**
	 * Joins a collection of strings into one string
	 * 
	 * @param all		the collection of substrings
	 * @param joinby	the token that joins the substrings
	 * @return			the final product: str1 + joinby + str2 + . . . + strN
	 */	
	public static String StringJoin(double[] all, String joinby){
		if (all == null){
			return " NULL ARRAY ";
		}		
		String joined = "";
		for (int i = 0; i < all.length; i++){
			joined += all[i];
			if (i < all.length - 1)
				joined += joinby;
		}
		return joined;
	}
	
	/**
	 * Joins a collection of strings into one string
	 * 
	 * @param all		the collection of substrings
	 * @param joinby	the token that joins the substrings
	 * @return			the final product: str1 + joinby + str2 + . . . + strN
	 */	
	public static String StringJoin(int[] all, String joinby){
		if (all == null){
			return " NULL ARRAY ";
		}		
		String joined = "";
		for (int i = 0; i < all.length; i++){
			joined += all[i];
			if (i < all.length - 1)
				joined += joinby;
		}
		return joined;
	}
	

	public static String StringJoin(BitSet all, String joinby) {
		if (all == null){
			return " NULL ARRAY ";
		}		
		String joined = "";
		int n = all.length();
//		System.out.println("n = " + n);
		for (int i = 0; i < n; i++){
			joined += all.get(i) ? "1" : "0";
			if (i < n - 1)
				joined += joinby;
		}
		return joined;
	}	
	
	/**
	 * Joins a collection of strings into one string with commas
	 * 
	 * @param all		the collection of substrings
	 * @return			the final product: str1 + joinby + str2 + . . . + strN
	 */	
	public static String StringJoin(boolean[] all){
		int[] all_ints = new int[all.length];
		for (int i = 0; i < all.length; i++){
			all_ints[i] = all[i] ? 1 : 0;
		}
		return StringJoin(all_ints, ", ");
	}
	
	/**
	 * Joins a collection of strings into one string with commas
	 * 
	 * @param all		the collection of substrings
	 * @return			the final product: str1 + joinby + str2 + . . . + strN
	 */	
	public static String StringJoin(int[] all){
		return StringJoin(all, ", ");
	}
	
	/**
	 * Joins a collection of strings into one string with commas
	 * 
	 * @param all		the collection of substrings
	 * @return			the final product: str1 + joinby + str2 + . . . + strN
	 */	
	public static String StringJoin(double[] all){
		return StringJoin(all, ", ");
	}

	/**
	 * Joins a collection of strings into one string with commas
	 * 
	 * @param all		the collection of substrings
	 * @return			the final product: str1 + joinby + str2 + . . . + strN
	 */	
	public static String StringJoin(ArrayList<Object> all){
		return StringJoin(all, ", ");
	}	
	
	/**
	 * Joins a collection of strings into one string
	 * 
	 * @param all		the collection of substrings
	 * @param joinby	the token that joins the substrings
	 * @return			the final product: str1 + joinby + str2 + . . . + strN
	 */	
	public static String StringJoin(Object[] all, String joinby){
		String joined = "";
		for (int i = 0; i < all.length; i++){
			joined += all[i];
			if (i < all.length - 1)
				joined += joinby;
		}
		return joined;
	}	
	
	/**
	 * Joins a collection of strings into one string with commas
	 * 
	 * @param all		the collection of substrings
	 * @return			the final product: str1 + joinby + str2 + . . . + strN
	 */	
	public static String StringJoin(Object[] all){
		return StringJoin(all, ", ");
	}	
	
	/**
	 * Joins a collection of strings into one string
	 * 
	 * @param all		the collection of substrings
	 * @param joinby	the token that joins the substrings
	 * @return			the final product: str1 + joinby + str2 + . . . + strN
	 */	
	public static String StringJoinStrings(Collection<String> all, String joinby){
		Object[] arr = all.toArray();
		String joined = "";
		for (int i = 0; i < arr.length; i++){
			joined += (String)arr[i];
			if (i < arr.length - 1)
				joined += joinby;
		}
		return joined;
	}
	
	/**
	 * Joins a collection of strings into one string with commas
	 * 
	 * @param all		the collection of substrings
	 * @return			the final product: str1 + joinby + str2 + . . . + strN
	 */	
	public static String StringJoin(Collection<String> all){
		return StringJoinStrings(all, ", ");
	}	
	
	/**
	 * Returns the max of a vector
	 * 
	 * @param values	The values of interest
	 * @return			The maximum of those values
	 */
    public static double max(double[] values) {
    	double max = Double.NEGATIVE_INFINITY;
        for (double value : values) {
        	if (value > max){
        		max = value;
        	}
        }
        return max;
    }
    
    /**
     * Sums an array
     * 
     * @param arr	The values of interest
     * @return		The sum of those values
     */
    public static double sum_array(double[] arr){
    	double sum = 0;
    	for (int i = 0; i < arr.length; i++){
    		sum += arr[i];
    	}
    	return sum;
    }
    /**
     * Sums an array
     * 
     * @param arr	The values of interest
     * @return		The sum of those values
     */
    public static int sum_array(int[] arr){
    	int sum = 0;
    	for (int i = 0; i < arr.length; i++){
    		sum += arr[i];
    	}
    	return sum;
    }    
    
    /**
     * Sums the inverse values of an array
     * 
     * @param arr	The values of interest
     * @return		The sum of the inverses of those values
     */
	public static double sum_inv_array(double[] arr) {
    	double sum = 0;
    	for (int i = 0; i < arr.length; i++){
    		sum += 1 / arr[i];
    	}
    	return sum;
	}	    
 
	/**
	 * Normalizes an array by dividing each value by the array's sum
	 * 
	 * @param arr	The values of interest
	 */
    public static void normalize_array(double[] arr){
    	double weight = sum_array(arr);
    	for (int i = 0; i < arr.length; i++){
    		arr[i] = arr[i] / weight;
    	}
    }
    	
	/**
	 * Weights an array by dividing each value by a specified value
	 * 
	 * @param weight	The value to divide each value in the array by
	 * @param arr		The values of interest
	 */
    public static void weight_arr(double[] arr, double weight){
    	for (int i = 0; i < arr.length; i++){
    		arr[i] = arr[i] / weight;
    	}
    }    

    /**
     * Subtracts one array from the other
     * 
     * @param arr1	The array of minuends
     * @param arr2	The array of subtrahends
     * @return		The array of differences
     */
	public static double[] subtract_arrays(double[] arr1, double[] arr2) {
		int n = arr1.length;
		double[] diff = new double[n];
		for (int i = 0; i < n; i++){
			diff[i] = arr1[i] - arr2[i];
		}
		return diff;
	}

    /**
     * Adds one array to another
     * 
     * @param arr1	The array of first addends
     * @param arr2	The array of seconds addends
     * @return		The array of sums
     */
	public static double[] add_arrays(double[] arr1, double[] arr2) {
		int n = arr1.length;
		double[] sum = new double[n];
		for (int i = 0; i < n; i++){
			sum[i] = arr1[i] + arr2[i];
		}
		return sum;
	}

    /**
     * Adds one array to another
     * 
     * @param arr1	The array of first addends
     * @param arr2	The array of seconds addends
     * @return		The array of sums
     */	
	public static double[] add_arrays(double[] arr1, int[] arr2) {
		int n = arr1.length;
		double[] sum = new double[n];
		for (int i = 0; i < n; i++){
			sum[i] = arr1[i] + arr2[i];
		}
		return sum;
	}

    /**
     * Adds one array to another
     * 
     * @param arr1	The array of first addends
     * @param arr2	The array of seconds addends
     * @return		The array of sums
     */		
	public static int[] add_arrays(int[] arr1, int[] arr2) {
		int n = arr1.length;
		int[] sum = new int[n];
		for (int i = 0; i < n; i++){
			sum[i] = arr1[i] + arr2[i];
		}
		return sum;
	}

    /**
     * Adds one array to another after first converting each addend to binary
     * (1 if the value > 0, 0 otherwise)
     * 
     * @param arr1	The array of first addends
     * @param arr2	The array of seconds addends
     * @return		The array of sums of binary valus
     */
	public static int[] binary_add_arrays(int[] arr1, int[] arr2) {
		int n = arr1.length;
		int[] sum = new int[n];
		for (int i = 0; i < n; i++){
			sum[i] = (arr1[i] >= 1 ? 1 : 0) + (arr2[i] >= 1 ? 1 : 0);
		}
		return sum;
	}

	public static int min_index(double[] vals) {
		int index_of_min = -99;
		double max_so_far = Double.MAX_VALUE;
		for (int i = 0; i < vals.length; i++){
			if (vals[i] < max_so_far){
				max_so_far = vals[i];
				index_of_min = i;
			}
		}
		return index_of_min;
	}

	public static int[] convert_bitvector_to_intvector(BitSet bitvector, int n) {
		int[] intvector = new int[n];
		for (int i = 0; i < n; i++){
			intvector[i] = bitvector.get(i) ? 1 : 0;
		}
		return intvector;
	}

	public static int[] convertArrayListToVecInt(ArrayList<Integer> arr) {
		int[] vec = new int[arr.size()];
		for (int i = 0; i < vec.length; i++){
			vec[i] = arr.get(i);
		}
		return vec;
	}

	public static int[] findIndicies(int[] vec, int val) {
		ArrayList<Integer> indicies = new ArrayList<Integer>();
		for (int i = 0; i < vec.length; i++){			
			if (vec[i] == val){				
				indicies.add(i);
			}
		}
		return convertArrayListToVecInt(indicies);
	}

	public static int[] convertIntegerListToPrimVec(Integer[] arr) {
		int[] vec = new int[arr.length];
		for (int i = 0; i < arr.length; i++){
			vec[i] = arr[i];
		}
		return vec;
	}	

}
