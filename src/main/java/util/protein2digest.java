package main.java.util;


import com.compomics.util.experiment.biology.enzymes.Enzyme;
import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.identification.protein_sequences.digestion.ExtendedPeptide;
import com.compomics.util.experiment.identification.protein_sequences.digestion.IteratorFactory;
import com.compomics.util.experiment.identification.protein_sequences.digestion.SequenceIterator;
import com.compomics.util.parameters.identification.search.DigestionParameters;
import main.java.pg.DatabaseInput;
import net.sf.jfasta.FASTAElement;
import net.sf.jfasta.FASTAFileReader;
import net.sf.jfasta.impl.FASTAElementIterator;
import net.sf.jfasta.impl.FASTAFileReaderImpl;
import org.apache.commons.cli.*;
import org.apache.commons.lang3.StringUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;


public class protein2digest {


    public static void main(String[] args) throws IOException, ParseException, InterruptedException {

        Cloger.getInstance().logger.info("Start analysis");
        Cloger.getInstance().logger.info(StringUtils.join(args," "));

        Options options = new Options();

        options.addOption("db", true, "Fasta format database file");
        options.addOption("e",true,"Enzyme, default is trypsin");
        options.addOption("minLength",true,"Min length of peptide");
        options.addOption("maxLength",true,"Max length of peptide");
        options.addOption("minMass",true,"Min mass of peptide");
        options.addOption("maxMass",true,"Max mass of peptide");
        options.addOption("miss",true,"The number of miss cleavage, default is 2");
        options.addOption("o",true,"Output file");

        options.addOption("h", false, "Help");

        options.addOption("v", false, "Debug");

        CommandLineParser parser = new PosixParser();
        CommandLine cmd = parser.parse(options, args);
        if (cmd.hasOption("h") || cmd.hasOption("help") || args.length == 0) {
            HelpFormatter f = new HelpFormatter();

            f.printHelp("Options", options);
            System.exit(0);
        }

        String db = null;
        if(cmd.hasOption("db")){
            db = cmd.getOptionValue("db");
        }

        int miss = 2;
        if(cmd.hasOption("miss")){
            miss = Integer.valueOf(cmd.getOptionValue("miss"));
        }

        String outfile = null;
        if(cmd.hasOption("o")){
            outfile = cmd.getOptionValue("o");
        }

        protein2digest pro2digest = new protein2digest();

        if(cmd.hasOption("minLength")){
            pro2digest.minLength = Integer.valueOf(cmd.getOptionValue("minLength"));
        }

        if(cmd.hasOption("maxLength")){
            pro2digest.maxLength = Integer.valueOf(cmd.getOptionValue("maxLength"));
        }

        if(cmd.hasOption("maxMass")){
            pro2digest.maxMass = Double.valueOf(cmd.getOptionValue("maxMass"));
        }

        if(cmd.hasOption("minMass")){
            pro2digest.minMass = Double.valueOf(cmd.getOptionValue("minMass"));
        }

        pro2digest.digest(db,Integer.parseInt(cmd.getOptionValue("e")),miss,outfile);

    }


    public double minMass = 600.0;
    public double maxMass = 6000.0;
    public int minLength = 7;
    public int maxLength = 50;



    private void digest(String db, int enzyme_index, int nmiss, String outfile) throws IOException, InterruptedException {
        File dbFile = new File(db);
        FASTAFileReader reader = new FASTAFileReaderImpl(dbFile);
        FASTAElementIterator it = reader.getIterator();
        // digest protein
        ArrayList<String> fixedModifications = new ArrayList<>();
        Enzyme enzyme = DatabaseInput.getEnzymeByIndex(enzyme_index);
        DigestionParameters digestionParameters = DatabaseInput.getDigestionPreferences(enzyme.getName(), nmiss);
        IteratorFactory iteratorModifications = new IteratorFactory(fixedModifications);

        BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(new File(outfile)));

        bufferedWriter.write("acc\tpeptide\tposition\tmiss\tmass\n");

        int num = 0;
        while (it.hasNext()) {
            FASTAElement el = it.next();
            el.setLineLength(1);
            String headLine[] = el.getHeader().split("\\s+");
            String proID = headLine[0];
            num++;
            String proSeq = el.getSequence().toUpperCase();

            proSeq = proSeq.replaceAll("\\*", "");

            if((num % 1000)==0){
                System.out.println(num);
            }

            // digest
            SequenceIterator sequenceIterator = iteratorModifications.getSequenceIterator(proSeq, digestionParameters, this.minMass, this.maxMass);

            //for(ProteinSequenceIterator.PeptideWithPosition peptideWithPosition: peptides){
            ExtendedPeptide peptideWithPosition;
            while ((peptideWithPosition = sequenceIterator.getNextPeptide()) != null) {
                Peptide peptide = peptideWithPosition.peptide;
                if(peptide.getSequence().length() >= this.minLength && peptide.getSequence().length() <= this.maxLength){
                    bufferedWriter.write(proID + "\t" + peptide.getSequence() + "\t" + peptideWithPosition.position + "\t" + peptide.getNMissedCleavages(digestionParameters) + "\t" + peptide.getMass() + "\n");

                }

            }


        }
        reader.close();
        bufferedWriter.close();

        Cloger.getInstance().logger.info("Proteins in total:"+num);
    }


}
