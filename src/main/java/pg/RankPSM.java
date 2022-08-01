package main.java.pg;

import org.apache.commons.lang3.StringUtils;

import java.io.*;
import java.util.*;

/**
 * Sort PSM result table and add a column (rank)
 */
public class RankPSM {

    public double pvalue;
    public double score;
    public double tol;
    /**
     * A string for sort when both pvalue and score are the same for PSMs.
     */
    public String string_for_sort;

    public String line;


    public RankPSM(double p,double s, double tol, String string_for_sort){
        this.pvalue = p;
        this.score = s;
        this.tol = Math.abs(tol);
        this.string_for_sort = string_for_sort;
    }

    public RankPSM(){

    }


    public  static String rankPSMs(String input) throws IOException {
        BufferedReader bReader = new BufferedReader(new FileReader(new File(input)));

        String head = bReader.readLine();
        head = head.trim();
        String []h = head.split("\t");
        HashMap<String,Integer> hIndex = new HashMap<>();
        for(int i=0;i<h.length;i++){
            hIndex.put(h[i],i);
        }
        String line;

        HashMap<String,ArrayList<RankPSM>> psmMap = new HashMap<>();
        while((line = bReader.readLine())!=null){
            line = line.trim();
            String []d = line.split("\t");

            double pvalue = Double.parseDouble(d[hIndex.get("pvalue")]);
            double score = Double.parseDouble(d[hIndex.get("score")]);
            double tol = Double.parseDouble(d[hIndex.get("tol_ppm")]);
            String msID = d[hIndex.get("spectrum_title")];

            //RankPSM rankPSM = new RankPSM(pvalue,score,line);
            RankPSM rankPSM = new RankPSM(pvalue,score,tol,d[hIndex.get("peptide")]+"|"+d[hIndex.get("modification")]);
            rankPSM.line = line;
            if(psmMap.containsKey(msID)){
                psmMap.get(msID).add(rankPSM);
            }else{

                ArrayList<RankPSM> psms = new ArrayList<>();
                psms.add(rankPSM);
                psmMap.put(msID,psms);

            }
        }

        String out = input;
        if(out.endsWith(".txt")){
            out = out.replaceFirst(".txt$","_rank.txt");
        }else{
            out = out.concat("_rank.txt");
        }

        BufferedWriter bWriter = new BufferedWriter(new FileWriter(new File(out)));
        bWriter.write(head+"\trank\n");
        if(psmMap.size()>=1){
            for(String msID : psmMap.keySet()){
                ArrayList<RankPSM> psms = psmMap.get(msID);
                RankPSMCompare vc = new RankPSMCompare();
                Collections.sort(psms,vc);
                for(int i=0;i<psms.size();i++){
                    int rank = i+1;
                    String []d = psms.get(i).line.split("\t");
                    if(rank>=2) {
                        d[hIndex.get("n_db")] = "-1";
                        d[hIndex.get("total_db")] = "-1";
                    }
                    bWriter.write(StringUtils.join(d,"\t") +"\t"+rank+"\n");
                }
            }
        }

        bWriter.close();
        bReader.close();
        return out;

    }

    /**
     * sort RankPSM: pvalue and score
     */
    public static class RankPSMCompare implements Comparator<RankPSM> {


        public int compare(RankPSM o1, RankPSM o2) {
            double p1 = o1.pvalue;
            double p2 = o2.pvalue;
            double s1 = o1.score;
            double s2 = o2.score;
            double t1 = o1.tol;
            double t2 = o2.tol;
            String c1 = o1.string_for_sort;
            String c2 = o2.string_for_sort;
            double r1 = p2 - p1;
            double r2 = s2 - s1;
            double r3 = t2 - t1;
            int r4 = c1.compareTo(c2);
            // score
            if(r2 > 0){
                return 1;
            }else if(r2 == 0){
                // p-value
                if(r1 > 0){
                    return -1;
                }else if(r1==0){
                    if(r3 > 0){
                        return -1;
                    }else if(r3==0){
                        if(r4 > 0){
                            return 1;
                        }else{
                            return -1;
                        }
                    }else{
                        return 1;
                    }
                }else{
                    return 1;
                }

            }else{
                return -1;
            }
        }

    }

    public static void main(String[] args) throws IOException {

        RankPSM.rankPSMs(args[0]);
    }


}
