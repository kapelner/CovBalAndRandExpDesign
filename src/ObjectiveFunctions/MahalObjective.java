package ObjectiveFunctions;

public class MahalObjective extends ObjectiveFunction {

	private double[][] Sinvmat;
	private double prop_const;

	public MahalObjective(double[][] Sinvmat, int n) {
		this.Sinvmat = Sinvmat;
		this.prop_const =  n / (double)4; // n * p_w * (1 - p_w) where p_w is #1/n
	}

	@Override
	public double calc(boolean debug_mode) {
		//as.numeric(t(X_T_bar_minus_X_C_bar) %*% inv_cov_X %*% X_T_bar_minus_X_C_bar)
		int p = XTbar.length;
//		System.out.println("p = " + p);
		double[] X_T_bar_minus_X_C_bar = new double[p];
//		System.out.println("X_T_bar_minus_X_C_bar.toString() = " + X_T_bar_minus_X_C_bar.toString());
		for (int j = 0; j < p; j++){
			X_T_bar_minus_X_C_bar[j] = XTbar[j] - XCbar[j];
		}
		
		double[] Sinvmat_times_X_T_bar_minus_X_C_bar = new double[p];
		
		//matrix multiplied by vector because a vector
		for (int j = 0; j < p; j++){
			double dot_product = 0;
			//the inner loop is just a vector times a vector which is a dot product
			for (int i = 0; i < p; i++){
				dot_product += Sinvmat[j][i] * X_T_bar_minus_X_C_bar[i]; //row fixed iterate over column
			}	
			Sinvmat_times_X_T_bar_minus_X_C_bar[j] = dot_product;
		}
		
		//vector times a vector is a dot product
		double dot_product = 0;
		for (int j = 0; j < p; j++){
			dot_product += Sinvmat_times_X_T_bar_minus_X_C_bar[j] * X_T_bar_minus_X_C_bar[j];
		}
		return prop_const * dot_product;
	}

}
