package ExperimentalDesign;

import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import CustomLogging.FileLoggedClass;
import ObjectiveFunctions.ObjectiveFunction;
import gurobi.GRBException;

public abstract class AllExperimentalDesigns extends FileLoggedClass {
	
	public Random rand_obj;
	
	//set by user
	protected int n;
	protected int p;
	protected String objective;
	protected Integer num_cores;
	protected boolean search_stopped;
	
	//data inputed from the user's data
	protected double[][] X;
	protected double[][] Sinv;	
	protected boolean wait;
	
	//temporary objects needed for search
	protected ExecutorService search_thread_pool;	
	protected boolean began_search;
	protected long t0;
	protected Long tf;
	
	
	public AllExperimentalDesigns(){
		num_cores = 1;
		rand_obj = new Random();
	}	
	
	public void beginSearch() {
//		System.out.println("beginSearch");
		began_search = true;
		
		t0 = System.currentTimeMillis();
		//build the pool and all tasks to it
		search_thread_pool = Executors.newFixedThreadPool(num_cores == null ? 1 : num_cores);
	}
	
	
	protected void afterBeginSearch() {
		search_thread_pool.shutdown(); //run em all (but not on this thread!)
		Thread await_completion = new Thread(){
			public void run(){
				try {
					search_thread_pool.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS); //infinity
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
				tf = System.currentTimeMillis();
			}
		};
		await_completion.start();
		if (wait){
			try {
				await_completion.join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}	
	

	public int timeElapsedInSeconds(){
		if (tf == null){
			return (int)(System.currentTimeMillis() - t0) / 1000;
		}
		return (int)(tf - t0) / 1000;
	}
	
	public long timeFinished(){
		return tf;
	}
	
	public boolean began(){
		return began_search;
	}
	
	public void stopSearch(){
		search_stopped = true;
	}	
	
	public void setObjective(String objective) throws Exception{
		if (!ObjectiveFunction.isValidObjFunction(objective)){
			throw new Exception("objective function not recognized");
		}
		this.objective = objective;
	}
	
	public void setNumCores(int num_cores){
//		System.out.println("setNumCores " +num_cores);
		this.num_cores = num_cores;
	}	
	
	public void setNandP(int n, int p) throws Exception{
		if (n % 2 != 0){
			throw new Exception("n must be even");
		}
//		System.out.println("setNandP n = " + n + " p = " + p);
		this.n = n;
		this.p = p;
	}	
	
	public void setDataRow(int i0, double[] x_i){
//		System.out.println("setDataRow " + i0 + "  " + x_i);
		if (X == null){
			X = new double[n][p];
		}
		for (int j = 0; j < p; j++){
			X[i0][j] = x_i[j];
		}
	}
	
	public void setDataRow(int i0, double x_i){
//		System.out.println("setDataRow " + i0 + "  " + x_i);
		if (X == null){
			X = new double[n][p];
		}
		double[] row = {x_i};
		X[i0] = row;
	}	
	
	public void setInvVarCovRow(int j0, double[] Sinv_i){
//		System.out.println("setInvVarCovRow " + j0 + "  " + Sinv_i);
		if (Sinv == null){
			Sinv = new double[p][p];
		}
		for (int j = 0; j < p; j++){
			Sinv[j0][j] = Sinv_i[j];
		}
	}
	
	public void setInvVarCovRow(int j0, double Sinv_i){
//		System.out.println("setInvVarCovRow " + j0 + "  " + Sinv_i);
		if (Sinv == null){
			Sinv = new double[p][p];
		}
		double[] row = {Sinv_i};
		Sinv[j0] = row;
	}
	
	public void setSeed(int seed){
		rand_obj.setSeed(seed);
	}	
	
	public void setWait(){
		wait = true;
	}	
}
