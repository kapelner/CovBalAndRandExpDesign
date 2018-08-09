package KarpExperimentalDesign;

import java.util.Collections;
import java.util.Comparator;


public class KarpDesignSearcherBalanced extends KarpDesignSearcher {
	
	public KarpDesignSearcherBalanced(double[][] Xstd) {
		super(Xstd);
	}
	
	
	protected class ObsBundleCompareBalanced implements Comparator<ObsBundle> {
		@Override
		public int compare(ObsBundle o1, ObsBundle o2) {
			//we cannot merge a singleton with anything else, so handle that here
			if (o1.size() == 1 || o2.size() == 1){
				if (o1.size() < o2.size()){
					return -1;
				}
				else if (o1.size() > o2.size()){
					return 1;
				}
				//if equal... then move to next return statement
			}
			//otherwise, we just compare based on value
			return Double.compare(o2.x_val, o1.x_val);
		}	
	}
	
	public void sortObsBundles(){
		Collections.sort(obs_bundles, new ObsBundleCompareBalanced());
	}
}
