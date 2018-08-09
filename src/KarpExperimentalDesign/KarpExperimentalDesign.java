package KarpExperimentalDesign;

import ExperimentalDesign.AllExperimentalDesigns;
import ExperimentalDesign.Tools;
import ObjectiveFunctions.*;

public class KarpExperimentalDesign extends AllExperimentalDesigns {
	
	
	private boolean balanced;
	private KarpDesignSearcher keds;
	
	//running the Java as standalone is for debug purposes ONLY!!!
	public static void main(String[] args) throws Exception{	

		for (int g = 0; g < 1; g++){
			KarpExperimentalDesign kd = new KarpExperimentalDesign();
			//set seed here for reproducibility during debugging
			kd.rand_obj.setSeed(1984);
	
			
			
			double[] temp = {7,8,5,10,15,16,9,3,21,1,37,25,28,19,42,65,13,43,58,49,89,30,40,50};//{1, 3, 5, 7, 10, 20, 30, 35, 40, 45, 50, 55, 80, 100};
			int n = temp.length;
//			int n = 14;
			kd.setNandP(n, 1);
			for (int i = 0; i < n; i++){
				double[] x_i = new double[1];
				x_i[0] = kd.rand_obj.nextDouble();
				x_i[0] = n - (i);
				x_i[0] = temp[i];
				kd.setDataRow(i, x_i);
			}
	//		System.out.println("Xstd");
	//		for (int i = 0; i < n; i++){
	//			System.out.println(Tools.StringJoin(od.Xstd[i]));
	//		}
			kd.setObjective(ObjectiveFunction.ABS);
			kd.setWait();
			kd.setBalanced();
			kd.beginSearch();
		}
//		System.out.println("progress: " + od.progress());
	}	
	
	public void beginSearch(){
		super.beginSearch();
    	search_thread_pool.execute(new Runnable(){
			public void run() {
				keds = balanced ? new KarpDesignSearcherBalanced(X) : new KarpDesignSearcherUnbalanced(X);
			}
		});
		afterBeginSearch();	
		System.out.println("FINAL INDIC_T: " + Tools.StringJoin(getKarpIndicT()));
		System.out.println("Num T: " + Tools.sum_array(getKarpIndicT()) + " n: " + n);
		System.out.println("Final obj val: " + getKarpObjectiveVal());
	}
	
	public void setBalanced(){
		balanced = true;
	}
	
	public double getKarpObjectiveVal(){		
		return keds.getObjVal();
	}
	
	public int[] getKarpIndicT() {
		return keds.getIndicT();
	}	
	
	public double progress(){
		return keds == null ? 0 : keds.progress();
	}	
}
