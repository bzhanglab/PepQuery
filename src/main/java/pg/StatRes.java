package main.java.pg;


import java.io.IOException;

import static main.java.pg.PeptideSearchMT.summariseResult;

/**
 * Get summary information (psm_rank.txt, ptm.txt)
 */
public class StatRes {

    public static void main(String[] args) throws IOException {
        if(args.length==3){
            CParameter.scoreMethod = Integer.valueOf(args[2]);
            summariseResult(args[0],args[1]);
        }else{
            System.err.println("<psm_rank.txt> <ptm.txt> <method>");
            System.exit(0);
        }

    }
}
