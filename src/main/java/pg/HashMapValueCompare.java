package main.java.pg;

import java.util.Comparator;
import java.util.HashMap;

/**
 * Sort a hashmap based on value => low to high
 *
 */

public class HashMapValueCompare implements Comparator<HashMap.Entry<String, Integer>> {


    public int compare(HashMap.Entry<String, Integer> o1, HashMap.Entry<String, Integer> o2) {
        int v1 = o1.getValue();
        int v2 = o2.getValue();
        int result = v2 - v1;
        if(result > 0){
            return -1;
        }else if(result==0){
            return 0;
        }else{
            return 1;
        }
    }

}
