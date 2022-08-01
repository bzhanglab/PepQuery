package main.java.pg;

import main.java.download.FileDownload;
import main.java.index.BuildMSlibrary;
import main.java.raw.RawDataConvert;
import main.java.util.protein2digest;

public class IMain {

    public static void main(String[] args) throws Exception {

        if(args.length<=0){
            PeptideSearchMT.main(args);
        }else{
            String program = args[0];
            String parameters[] = new String[args.length-1];
            for(int i=0;i<parameters.length;i++){
                parameters[i] = args[i+1];
            }

            if(program.equalsIgnoreCase("digest")){
                protein2digest.main(parameters);
            }else if(program.equalsIgnoreCase("index")){
                BuildMSlibrary.main(parameters);
            }else if(program.equalsIgnoreCase("download")) {
                FileDownload.main(parameters);
            }else if(program.equalsIgnoreCase("MSConvert")){
                RawDataConvert.main(parameters);
            }else{
                PeptideSearchMT.main(args);
            }
        }


    }
}
