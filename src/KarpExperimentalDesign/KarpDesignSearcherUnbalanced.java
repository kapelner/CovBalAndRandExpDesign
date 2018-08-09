package KarpExperimentalDesign;

import java.util.Collections;
import java.util.Comparator;

public class KarpDesignSearcherUnbalanced extends KarpDesignSearcher {
	
	public KarpDesignSearcherUnbalanced(double[][] Xstd) {
		super(Xstd);
	}
	
	protected class ObsBundleCompareUnbalanced implements Comparator<ObsBundle> {
		@Override
		public int compare(ObsBundle o1, ObsBundle o2) {
			return Double.compare(o2.x_val, o1.x_val);
		}		
	}	
	
	public void sortObsBundles(){
		Collections.sort(obs_bundles, new ObsBundleCompareUnbalanced());
	}
}
