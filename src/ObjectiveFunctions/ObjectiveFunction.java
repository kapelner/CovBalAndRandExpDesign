package ObjectiveFunctions;

import java.util.ArrayList;

public abstract class ObjectiveFunction {

	//valid objective functions
	public static final String MAHAL = "mahal_dist";
	public static final String ABS = "abs_sum_diff";
	private static final ArrayList<String> VALID_OBJ_FUNCTIONS = new ArrayList<String>();
	static {
		VALID_OBJ_FUNCTIONS.add(MAHAL);
		VALID_OBJ_FUNCTIONS.add(ABS);
	};
	
	protected double[] XTbar;
	protected double[] XCbar;

	public abstract double calc(boolean debug_mode);
	
	public static boolean isValidObjFunction(String objective){
		return VALID_OBJ_FUNCTIONS.contains(objective);
	}

	public void setXTbar(double[] XTbar){
//		System.out.println("XTbar: " + Tools.StringJoin(XTbar.getData()));
		this.XTbar = XTbar;
	}
	
	public void setXCbar(double[] XCbar){
//		System.out.println("XCbar: " + Tools.StringJoin(XCbar.getData()));
		this.XCbar = XCbar;
	}
}
