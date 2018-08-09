package KarpExperimentalDesign;

import java.util.ArrayList;
import java.util.LinkedList;

import ObjectiveFunctions.*;
import ExperimentalDesign.Tools;

public abstract class KarpDesignSearcher {
	
	protected Integer[] indicT;
	protected int n;
	protected LinkedList<ObsBundle> obs_bundles;
	protected double[][] Xstd;
	//first allocation is always treatment
	protected static final int FIRST_ALLOCATION = 1;
	
	protected class ObsBundle {
		public ObsBundle home;
		public ObsBundle a;
		public ObsBundle b;		
		public Double x_val;
		public Integer indicT_index;
		
		//only used for the initial population of observation bundles
		public ObsBundle(Double x_val, Integer indicT_index) {
			this.x_val = x_val;
			this.indicT_index = indicT_index;
		}
		
		//only used for merge!
		public ObsBundle(ObsBundle a, ObsBundle b) {
			this.a = a;
			a.home = this;
			this.b = b;
			b.home = this;
			//now we gotta set the values for the first one
			if (a.simple() && a.allocationNotSet()){
				a.allocate(FIRST_ALLOCATION);
			}
			if (b.simple() && b.allocationNotSet()){
				b.allocate(1 - FIRST_ALLOCATION);
			}
			//b always gets flipped. a does not!
			else {
				b.flipAllocation();
			}

			x_val = updateXVal();
		}
		
		private void allocate(int alloc) {
			indicT[indicT_index] = alloc;
		}
		
		private boolean simple(){
			return (a == null && b == null) ? true : false;
		}
		
		private boolean allocationNotSet(){
			return indicT[indicT_index] == null ? true : false;
		}
		
		public void flipAllocation(){
			if (simple()){
				if (allocationNotSet()){
					allocate(FIRST_ALLOCATION);
				}
				else {
					allocate(1 - indicT[indicT_index]);
				}
			}
			else {
				ObsBundle temp = a;
				a = b;
				b = temp;
			}
			if (a != null){
				a.flipAllocation();
			}
			if (b != null){
				b.flipAllocation();
			}
		}
		
		private double updateXVal() {
			if (!simple()){
				if (a.simple() && b.simple()){			
					x_val = indicT[a.indicT_index] == 1 ? (a.updateXVal() - b.updateXVal()) : (b.updateXVal() - a.updateXVal());
				}
				else {
					x_val = a.updateXVal() + b.updateXVal();
				}				
			}
			return x_val;
		}

		public int size(){
			if (a == null && b == null){
				return 1;
			} 
			else if (a == null){
				return 1 + b.size();
			}
			else if (b == null){
				return 1 + a.size();
			}
			else {
				return a.size() + b.size();
			}
		}
		
		//debug only
		public void print(String indent){
			System.out.print(indent + "x_val = " + x_val);
			if (indicT_index != null && indicT[indicT_index] != null){
				System.out.print("  [" +  (indicT[indicT_index] == 1 ? "Treatment]" : "Control]") + "\n");
				
			}
			else {
				System.out.print("\n");
			}
			if (a != null){
				a.print(indent + "  (a) ");
			}			
			if (b != null){
				b.print(indent + "  (b) ");
			}
		}
	}
	
	public KarpDesignSearcher(double[][] Xstd){
		this.Xstd = Xstd;
		
		n = Xstd.length;
		obs_bundles = new LinkedList<ObsBundle>();
		for (int i = 0; i < n; i++){
			obs_bundles.add(new ObsBundle(Xstd[i][0], i));
		}
		//initialize the integer vector --- all places are null first
		indicT = new Integer[n];
			
		int iter = 1;
		while (true){
			//the first thing to do is order these things up
			sortObsBundles();			
			
			System.out.println("\n\niter " + iter + " size of obs_bundles: " + obs_bundles.size() + "  ===========================================================================================================");
			for (int i = 0; i < obs_bundles.size(); i++){
				System.out.println("  BUNDLE #" + (i + 1));
				obs_bundles.get(i).print("    ");
			}			
			
			
			//now that we've ordered them, we can merge the first two - they will become a pair
			mergeFirstTwoObsBundles();
			
			//now check if there is only one thing left
			if (obs_bundles.size() == 1){
				break;
			}
			iter++;
		}
		System.out.println("\n\niter FINAL size of obs_bundles: " + obs_bundles.size() + "  ===========================================================================================================");
		for (int i = 0; i < obs_bundles.size(); i++){
			System.out.println("  BUNDLE #" + (i + 1));
			obs_bundles.get(i).print("    ");
		}	
	}
	
	protected abstract void sortObsBundles();

	protected void mergeFirstTwoObsBundles() {
		//remove the first two
		ObsBundle a = obs_bundles.removeFirst();
		ObsBundle b = obs_bundles.removeFirst();
		//create a merged bundle in the proper order. The constructor handles the merging (kinda hacky, but whatever)
		ObsBundle ab_or_ba = a.x_val >= b.x_val ? new ObsBundle(a, b) : new ObsBundle(b, a);
		//put that in its place
		obs_bundles.addFirst(ab_or_ba);
	}

	public double getObjVal() {		
		int[] indicT = getIndicT();
		int[] i_Ts = Tools.findIndicies(indicT, 1);
//		System.out.println("i_Ts " + Tools.StringJoin(i_Ts));
		int[] i_Cs = Tools.findIndicies(indicT, 0);
		ArrayList<double[]> XT = Tools.subsetMatrix(Xstd, i_Ts); 
		ArrayList<double[]> XC = Tools.subsetMatrix(Xstd, i_Cs);
		double[] avg_Ts = Tools.colAvg(XT, 1);
		double[] avg_Cs = Tools.colAvg(XC, 1);	
		AbsSumObjective obj_fun = new AbsSumObjective();
		obj_fun.setXTbar(avg_Ts);
		obj_fun.setXCbar(avg_Cs);
		return obj_fun.calc(false);
	}
	
	public int[] getIndicT() {
		return Tools.convertIntegerListToPrimVec(this.indicT);
	}
	
	public double progress() {
		return (n - obs_bundles.size() + 1) / (double)n;
	}	
}
