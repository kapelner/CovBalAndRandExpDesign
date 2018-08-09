package ObjectiveFunctions;

import ExperimentalDesign.Tools;

public class AbsSumObjectiveWithDiagnostics extends ObjectiveFunction {

	private double[] xbardiffjs;

	public double[] getXbardiffjs() {
		return xbardiffjs;
	}

	@Override
	public double calc(boolean debug_mode) {
		if (debug_mode){
			System.out.println("XTbar: " + Tools.StringJoin(XTbar, ","));
			System.out.println("XCbar: " + Tools.StringJoin(XCbar, ","));			
		}
		int p = XTbar.length;
		xbardiffjs = new double[p];
		double abs_sum = 0;
		for (int j = 0; j < p; j++){
			double diff = XTbar[j] - XCbar[j];
			xbardiffjs[j] = diff;
			abs_sum += (diff < 0.0 ? -diff : diff); //faster than Math.abs accd to JProfiler
		}
		return abs_sum;
	}

}
