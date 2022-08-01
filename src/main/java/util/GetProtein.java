package main.java.util;

import net.sf.jfasta.FASTAElement;
import net.sf.jfasta.FASTAFileReader;
import net.sf.jfasta.impl.FASTAElementIterator;
import net.sf.jfasta.impl.FASTAFileReaderImpl;
import org.apache.commons.cli.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Fetch proteins according to prefix of the protein IDs, generate two new files, one is the file containing the
 * interested proteins and the other one is the file containing the remaining proteins.
 */
public class GetProtein {


    public static void main(String[] args) throws IOException, ParseException {


        Options options = new Options();

        //options.addOption("o", true, "Output dir");
        options.addOption("db", true, "Fasta format database file");
        options.addOption("c",true,"Matched string");
        options.addOption("pos",true,"The position of the string in protein ID, 1=>begin;2=>end;3=>any");
        options.addOption("o",true,"The output folder");
        options.addOption("r",false,"regex");

        options.addOption("h", false, "Help");

        options.addOption("v", false, "Debug");

        CommandLineParser parser = new PosixParser();
        CommandLine cmd = parser.parse(options, args);
        if (cmd.hasOption("h") || cmd.hasOption("help") || args.length == 0) {
            HelpFormatter f = new HelpFormatter();

            // System.out.println(version);
            f.printHelp("Options", options);
            System.exit(0);
        }


        String db = cmd.getOptionValue("db");
        String prefix = cmd.getOptionValue("c");
        int pos = 1;
        if(cmd.hasOption("pos")){
            pos = Integer.valueOf(cmd.getOptionValue("pos"));
        }

        String outdir = "./";
        if(cmd.hasOption("o")){
            outdir = cmd.getOptionValue("o");
        }

        String outdb1 = outdir+"/target_protein.fasta";
        String outdb2 = outdir+"/remain_protein.fasta";

        BufferedWriter bufferedWriter1 = new BufferedWriter(new FileWriter(new File(outdb1)));
        BufferedWriter bufferedWriter2 = new BufferedWriter(new FileWriter(new File(outdb2)));

        int n1 = 0;
        int n2 = 0;

        File dbFile = new File(db);
        FASTAFileReader reader = new FASTAFileReaderImpl(dbFile);
        FASTAElementIterator it = reader.getIterator();
        while (it.hasNext()) {
            FASTAElement el = it.next();
            el.setLineLength(1);
            String headLine[] = el.getHeader().split("\\s+");
            String proID = headLine[0];
            String proSeq = el.getSequence().toUpperCase();

            if(cmd.hasOption("r")){
                String ps = cmd.getOptionValue("c");
                ps = "^" + ps + "$";

                Pattern pattern = Pattern.compile(ps);
                Matcher matcher = pattern.matcher(proID);

                if (matcher.find()) {
                    bufferedWriter1.write(">" + el.getHeader() + "\n");
                    bufferedWriter1.write(proSeq + "\n");
                    n1++;
                } else {
                    bufferedWriter2.write(">" + el.getHeader() + "\n");
                    bufferedWriter2.write(proSeq + "\n");
                    n2++;
                }

            }else if(pos==1) {

                if (proID.startsWith(prefix)) {
                    bufferedWriter1.write(">" + el.getHeader() + "\n");
                    bufferedWriter1.write(proSeq + "\n");
                    n1++;
                } else {
                    bufferedWriter2.write(">" + el.getHeader() + "\n");
                    bufferedWriter2.write(proSeq + "\n");
                    n2++;
                }
            }else if(pos==2){
                if (proID.endsWith(prefix)) {
                    bufferedWriter1.write(">" + el.getHeader() + "\n");
                    bufferedWriter1.write(proSeq + "\n");
                    n1++;
                } else {
                    bufferedWriter2.write(">" + el.getHeader() + "\n");
                    bufferedWriter2.write(proSeq + "\n");
                    n2++;
                }
            }else if(pos==3){
                if (proID.contains(prefix)) {
                    bufferedWriter1.write(">" + el.getHeader() + "\n");
                    bufferedWriter1.write(proSeq + "\n");
                    n1++;
                } else {
                    bufferedWriter2.write(">" + el.getHeader() + "\n");
                    bufferedWriter2.write(proSeq + "\n");
                    n2++;
                }
            }

        }
        reader.close();

        bufferedWriter1.close();
        bufferedWriter2.close();

        Cloger.getInstance().logger.info(n1+"\t"+n2);


    }

}
