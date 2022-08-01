package main.java.util;


import org.apache.commons.lang3.StringUtils;

import java.io.*;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Generate tag file from pFind result
 */
public class generateTags {

    public static void main(String[] args) throws IOException {

        generateTagsFromOpen_pFind(args[0],"");

    }


    /**
     *
     * @param input_file *.qry.res
     * @param out_file
     */
    public static void generateTagsFromOpen_pFind(String input_file, String out_file) throws IOException {

        BufferedReader br = new BufferedReader(new FileReader(new File(input_file)));

        ArrayList<String> peps = new ArrayList<>();
        String current_spectrum_title = "";
        String line;
        int index_number = 0;
        Pattern pattern = Pattern.compile("^\\d");
        while((line = br.readLine()) != null){
            if(line.startsWith("S")){

                if(peps.size()>0){
                    Cloger.getInstance().logger.info(current_spectrum_title+"\t"+ StringUtils.join(peps,";"));
                }

                peps.clear();

                index_number=index_number+1;
                current_spectrum_title = br.readLine().trim();

            }
            Matcher m = pattern.matcher(line);

            if(m.find()){
                String d[] = line.split("\t");
                if(d.length >= 5){
                    peps.add(d[1]);
                }
            }
        }

        br.close();
    }
}
