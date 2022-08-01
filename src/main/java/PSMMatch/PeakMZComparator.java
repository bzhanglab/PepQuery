package main.java.PSMMatch;

import java.util.Comparator;

/**
 * Sort a ArrayList containing JPeak objects according to the peak mz, from the minimum to the maximum
 *
 */
public class PeakMZComparator implements Comparator<JPeak> {

	
	public int compare(JPeak p1, JPeak p2) {
		double  mz1 = p1.getMz();
		double  mz2 = p2.getMz();
		if(mz1>mz2){
			return 1;
		}else if(mz1==mz2){
			return 0;
		}else{
			return -1;
		}
	}

}
