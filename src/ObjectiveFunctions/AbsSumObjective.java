package ObjectiveFunctions;

import ExperimentalDesign.Tools;

public class AbsSumObjective extends ObjectiveFunction {

	@Override
	public double calc(boolean debug_mode) {
		if (debug_mode){
			System.out.println("XTbar: " + Tools.StringJoin(XTbar, ","));
			System.out.println("XCbar: " + Tools.StringJoin(XCbar, ","));			
		}
		
		double abs_sum = 0;
		for (int j = 0; j < XTbar.length; j++){
			double diff = XTbar[j] - XCbar[j];
			abs_sum += (diff < 0.0 ? -diff : diff); //faster than Math.abs accd to JProfiler
		}		
		return abs_sum;
	}

}
