package main.java.util;

import com.compomics.util.experiment.biology.enzymes.Enzyme;
import com.compomics.util.experiment.biology.modifications.Modification;
import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.experiment.identification.protein_sequences.digestion.ExtendedPeptide;
import com.compomics.util.experiment.identification.protein_sequences.digestion.IteratorFactory;
import com.compomics.util.experiment.identification.protein_sequences.digestion.SequenceIterator;
import com.compomics.util.parameters.identification.search.DigestionParameters;
import main.java.PSMMatch.JPeptide;
import main.java.pg.*;
import net.sf.jfasta.FASTAElement;
import net.sf.jfasta.FASTAFileReader;
import net.sf.jfasta.impl.FASTAElementIterator;
import net.sf.jfasta.impl.FASTAFileReaderImpl;
import org.apache.commons.cli.*;
import org.apache.commons.lang3.StringUtils;

import org.sqlite.SQLiteConfig;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.sql.*;
import java.util.ArrayList;
import java.util.HashSet;


public class createProDB {


    public static void main(String[] args) throws ParseException, IOException, InterruptedException, SQLException, ClassNotFoundException {

        Cloger.getInstance().logger.info("Start analysis");
        Cloger.getInstance().logger.info(StringUtils.join(args," "));

        Options options = new Options();

        //options.addOption("o", true, "Output dir");
        options.addOption("db", true, "Fasta format database file");
        options.addOption("fixMod",true,"Fixed modification");
        options.addOption("varMod",true,"Variable modification");
        options.addOption("maxVar",true,"Max number of variable modifications, default is 3");
        options.addOption("printPTM",false,"Print PTMs");
        //options.addOption("miss",true,"Max number of missed cleavage sites, default is 2");

        options.addOption("e",true,"0:Non enzyme, 1:Trypsin (default), 2:Trypsin (no P rule), 3:Arg-C, 4:Arg-C (no P rule), 5:Arg-N, 6:Glu-C, 7:Lys-C");
        options.addOption("c",true,"The max missed cleavages, default is 2");

        options.addOption("minLength", true, "The minimum length of peptide to consider, default is 7");
        options.addOption("maxLength", true, "The maximum length of peptide to consider, default is 45");



        options.addOption("h", false, "Help");

        options.addOption("v", false, "Debug");

        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = parser.parse(options, args);
        if (cmd.hasOption("h") || cmd.hasOption("help") || args.length == 0) {
            HelpFormatter f = new HelpFormatter();

            // System.out.println(version);
            f.printHelp("Options", options);
            System.exit(0);
        }


        if(cmd.hasOption("printPTM")){
            ModificationGear.getInstance().printPTM();
            System.exit(0);
        }

        int maxVar = 3;
        if(cmd.hasOption("maxVar")){
            maxVar = Integer.valueOf(cmd.getOptionValue("maxVar"));
        }

        if(cmd.hasOption("e")){
            CParameter.enzyme = Integer.valueOf(cmd.getOptionValue("e"));
        }
        if(cmd.hasOption("c")){
            int missCleavages = Integer.valueOf(cmd.getOptionValue("c"));
            if(missCleavages<0){
                missCleavages = 0;
            }
            CParameter.maxMissedCleavages = missCleavages;
        }

        //int minPeptideLength = 7;
        if(cmd.hasOption("minLength")){
            CParameter.minPeptideLength = Integer.valueOf(cmd.getOptionValue("minLength"));
        }

        //int maxPeptideLength = 45;
        if(cmd.hasOption("maxLength")){
            CParameter.maxPeptideLength = Integer.valueOf(cmd.getOptionValue("maxLength"));
        }


        DatabaseInput databaseInput = new DatabaseInput("");
        if(cmd.hasOption("fixMod")){
            CParameter.fixMods = cmd.getOptionValue("fixMod");
            String ifixmod[] = cmd.getOptionValue("fixMod").split(",");
            for (String i : ifixmod) {
                int id = Integer.valueOf(i);
                if(id != 0) {
                    databaseInput.addFixedModifications(id);
                }
            }
        }
        if(cmd.hasOption("varMod")){
            CParameter.varMods = cmd.getOptionValue("varMod");
            String ivarmod[] = cmd.getOptionValue("varMod").split(",");
            for (String i : ivarmod) {
                int id = Integer.valueOf(i);
                databaseInput.addVarModifications(id);
            }
        }


        String db = cmd.getOptionValue("db");
        String sqldbfile = db+".sqldb";


        CParameter.print();

        File sql_file = new File(sqldbfile);
        if(sql_file.exists()){
            Cloger.getInstance().logger.info("File exists:"+sqldbfile);
            if(sql_file.delete()){
                Cloger.getInstance().logger.info("Delete file:"+sqldbfile);
            }
        }


        run(db,sqldbfile,databaseInput.varModifications,databaseInput.fixedModifications,maxVar,true);

        Cloger.getInstance().logger.info("End analysis");

    }

    /**
    public static void run(String db, String sqldbfile,
                           ArrayList<String> varModifications,
                           ArrayList<String> fixedModifications,
                           int maxVarMods,int nm) throws SQLException, IOException, ClassNotFoundException, InterruptedException {

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // SQL database
        SQLiteConfig config = new SQLiteConfig();
        config.setSynchronous(SQLiteConfig.SynchronousMode.OFF);

        Connection connection = DriverManager.getConnection("jdbc:sqlite:" + sqldbfile, config.toProperties());
        connection.setAutoCommit(false);
        Statement stmt = connection.createStatement();
        StringBuilder sqlBuilder = null;

        stmt.executeUpdate("drop table if exists prodb");

        // create table
        sqlBuilder = new StringBuilder();
        sqlBuilder.append("create table prodb ");
        sqlBuilder.append("(");
        sqlBuilder.append("peptide text").append(",");
        sqlBuilder.append("mass double").append(",");
        sqlBuilder.append("protein text").append(",");
        sqlBuilder.append("position integer").append(",");
        sqlBuilder.append("modification text").append(")");

        //System.out.println(sqlBuilder.toString());
        stmt.executeUpdate(sqlBuilder.toString());
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // read protein database
        File dbFile = new File(db);
        FASTAFileReader reader = new FASTAFileReaderImpl(dbFile);
        FASTAElementIterator it = reader.getIterator();

        // variable modifications
        ModificationFactory ptmFactory = ModificationFactory.getInstance();
        ArrayList<Modification> varPTMs = new ArrayList<>();
        for(String ptmName: varModifications){
            varPTMs.add(ptmFactory.getModification(ptmName));
        }

        // digest protein
        //DigestionPreferences digestionPreferences = DigestionPreferences.getDefaultPreferences();
        Enzyme enzyme = DatabaseInput.getEnzymeByIndex(CParameter.enzyme);
        DigestionParameters digestionPreferences = DatabaseInput.getDigestionPreferences(enzyme.getName(),CParameter.maxMissedCleavages);

        IteratorFactory iteratorModifications = new IteratorFactory(fixedModifications);


        // if a peptide is searched before, then it will be omitted.
        HashSet<String> searchedPeptides = new HashSet<>();

        int num = 0;
        while (it.hasNext()) {
            FASTAElement el = it.next();
            el.setLineLength(1);
            String headLine[] = el.getHeader().split("\\s+");
            String proID = headLine[0];
            num++;

            if( (num%1000)==0){
                System.out.println("Finished proteins:"+num);

            }

            String proSeq = el.getSequence().toUpperCase();

            proSeq = proSeq.replaceAll("\\*","");
            proSeq = proSeq.replaceAll("I","L");

            // digest
            SequenceIterator sequenceIterator = iteratorModifications.getSequenceIterator(proSeq, digestionPreferences, CParameter.minPeptideMass, CParameter.maxPeptideMass);

            //for(ProteinSequenceIterator.PeptideWithPosition peptideWithPosition: peptides){
            ExtendedPeptide peptideWithPosition;
            while ((peptideWithPosition = sequenceIterator.getNextPeptide()) != null) {

                Peptide peptide = peptideWithPosition.peptide;

                if(peptide.getSequence().length() < CParameter.minPeptideLength || peptide.getSequence().length() > CParameter.maxPeptideLength){
                    continue;
                }


                //System.out.println(peptide.getSequence()+"\t"+peptideWithPosition.getPosition());
                if(searchedPeptides.contains(peptide.getSequence())){
                    continue;
                }else{
                    searchedPeptides.add(peptide.getSequence());
                }
                int pos = peptideWithPosition.position;
                // generate peptide isoform (different ptm)
                ArrayList<Peptide> peptideIsoforms = PeptideInput.calcPeptideIsoforms(peptide,varPTMs,maxVarMods);
                peptideIsoforms.add(peptide);
                for(Peptide p: peptideIsoforms){
                    double peptideMass = JPeptide.getMass(p);//p.getMass();

                    ModificationMatch modificationMatches[] = p.getVariableModifications();
                    ArrayList<String> mods = new ArrayList<>();
                    if(modificationMatches != null) {
                        for (ModificationMatch modificationMatch : modificationMatches) {
                            mods.add(modificationMatch.getSite() + ":" + modificationMatch.getModification());
                        }
                    }

                    String modstring;
                    if(mods.size()>=1){
                        modstring = StringUtils.join(mods,";");
                    }else{
                        modstring = "-";
                    }

                    sqlBuilder = new StringBuilder();

                    sqlBuilder.append(" INSERT INTO prodb (peptide,mass,protein,position,modification) ");
                    sqlBuilder.append("VALUES (");
                    sqlBuilder.append("'").append(p.getSequence()).append("'").append(",");
                    sqlBuilder.append(peptideMass).append(",");
                    sqlBuilder.append("'").append(proID).append("'").append(",");
                    sqlBuilder.append(pos).append(",");
                    sqlBuilder.append("'").append(modstring).append("'").append(" )");
                    stmt.executeUpdate(sqlBuilder.toString());

                }

            }


        }

        reader.close();
        System.out.println("Protein sequences:"+num);
        System.out.println("Unique peptide sequences:"+searchedPeptides.size());

        connection.commit();
        stmt.executeUpdate("create index mass_index on prodb (mass)");
        connection.commit();

        stmt.close();
        connection.close();


    }**/

    public static void run(String db, String sqldbfile,
                           ArrayList<String> varModifications,
                           ArrayList<String> fixedModifications,
                           int maxVarMods,boolean writeOBJ) throws SQLException, IOException, ClassNotFoundException, InterruptedException {

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // SQL database
        SQLiteConfig config = new SQLiteConfig();
        config.setSynchronous(SQLiteConfig.SynchronousMode.OFF);

        Connection connection = DriverManager.getConnection("jdbc:sqlite:" + sqldbfile, config.toProperties());
        connection.setAutoCommit(false);
        Statement stmt = connection.createStatement();
        StringBuilder sqlBuilder = null;

        stmt.executeUpdate("drop table if exists prodb;");

        // create table
        sqlBuilder = new StringBuilder();
        sqlBuilder.append("create table prodb ");
        sqlBuilder.append("(");
        sqlBuilder.append("sequence text").append(",");
        sqlBuilder.append("mass double").append(",");
        sqlBuilder.append("peptide blob").append(",");
        sqlBuilder.append("protein text").append(",");
        sqlBuilder.append("position integer").append(",");
        sqlBuilder.append("modification text").append(",");
        sqlBuilder.append("pep_mod text").append(")");
        //System.out.println(sqlBuilder.toString());
        stmt.executeUpdate(sqlBuilder.toString());

        String psql = "INSERT INTO prodb (sequence,mass,peptide,protein,position,modification,pep_mod) VALUES(?,?,?,?,?,?,?);";

        PreparedStatement pstmt = connection.prepareStatement(psql);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // read protein database
        File dbFile = new File(db);
        FASTAFileReader reader = new FASTAFileReaderImpl(dbFile);
        FASTAElementIterator it = reader.getIterator();

        // variable modifications
        ModificationFactory ptmFactory = ModificationFactory.getInstance();
        ArrayList<Modification> varPTMs = new ArrayList<>();
        for(String ptmName: varModifications){
            varPTMs.add(ptmFactory.getModification(ptmName));
        }

        ArrayList<Modification> fixedPTMs = new ArrayList<>();
        for(String ptmName: fixedModifications){
            fixedPTMs.add(ptmFactory.getModification(ptmName));
        }

        // digest protein
        // DigestionPreferences digestionPreferences = DigestionPreferences.getDefaultPreferences();
        Enzyme enzyme = DatabaseInput.getEnzymeByIndex(CParameter.enzyme);
        DigestionParameters digestionPreferences = DatabaseInput.getDigestionPreferences(enzyme.getName(),CParameter.maxMissedCleavages);

        ArrayList<String> fixT = new ArrayList<>();
        IteratorFactory iteratorModifications = new IteratorFactory(fixT);


        // if a peptide is searched before, then it will be omitted.
        HashSet<String> searchedPeptides = new HashSet<>();

        int num = 0;
        double peptideMass = 0.0;
        while (it.hasNext()) {
            FASTAElement el = it.next();
            el.setLineLength(1);
            String headLine[] = el.getHeader().split("\\s+");
            String proID = headLine[0];


            if( (num%1000)==0){
                Cloger.getInstance().logger.info("Finished proteins:"+num+", "+searchedPeptides.size());

            }

            num++;


            // System.out.println(num+", " + searchedPeptides.size());


            String proSeq = el.getSequence().toUpperCase();

            proSeq = proSeq.replaceAll("\\*","");
            proSeq = proSeq.replaceAll("I","L");

            // digest
            SequenceIterator sequenceIterator = iteratorModifications.getSequenceIterator(proSeq, digestionPreferences, CParameter.minPeptideMass, CParameter.maxPeptideMass);

            //for(ProteinSequenceIterator.PeptideWithPosition peptideWithPosition: peptides){
            ExtendedPeptide peptideWithPosition;
            while ((peptideWithPosition = sequenceIterator.getNextPeptide()) != null) {

                Peptide peptide = peptideWithPosition.peptide;
                if(peptide.getSequence().length() < CParameter.minPeptideLength || peptide.getSequence().length() > CParameter.maxPeptideLength){
                    continue;
                }


                //System.out.println(peptide.getSequence()+"\t"+peptideWithPosition.getPosition());
                if(searchedPeptides.contains(peptide.getSequence())){
                    continue;
                }else{
                    searchedPeptides.add(peptide.getSequence());
                }
                int pos = peptideWithPosition.position;


                // Consider modification
                ArrayList<Peptide> peptideIsoforms = PeptideInput.calcPeptideIsoforms(peptide.getSequence());
                //peptideIsoforms.add(peptide);
                for(Peptide p: peptideIsoforms){
                    JPeptide jPeptide = new JPeptide(p);
                    peptideMass = jPeptide.getMass();
                    ModificationMatch modificationMatches[] = p.getVariableModifications();

                    ArrayList<String> mods = new ArrayList<>();
                    if(modificationMatches != null) {
                        for (ModificationMatch modificationMatch : modificationMatches) {
                            mods.add(modificationMatch.getSite() + ":" + modificationMatch.getModification());
                        }
                    }

                    String modstring;
                    if(mods.size()>=1){
                        modstring = StringUtils.join(mods,";");
                    }else{
                        modstring = "-";
                    }

                    ByteArrayOutputStream bos = new ByteArrayOutputStream();

                    ObjectOutputStream oos = new ObjectOutputStream(bos);

                    oos.writeObject(p);

                    //pstmt.setString(1,p.getSequenceWithLowerCasePtms());
                    pstmt.setString(1,p.getSequence());
                    pstmt.setDouble(2,peptideMass);
                    pstmt.setBytes(3,bos.toByteArray());
                    pstmt.setString(4,proID);
                    pstmt.setInt(5,pos);
                    pstmt.setString(6,modstring);
                    pstmt.setString(7,p.getSequence());
                    pstmt.executeUpdate();

                }

                /**
                // generate peptide isoform (different ptm)
                ArrayList<Peptide> peptideIsoforms = PeptideInput.calcPeptideIsoforms(peptide,varPTMs,maxVarMods);
                peptideIsoforms.add(peptide);
                for(Peptide p: peptideIsoforms){
                    double peptideMass = JPeptide.getMass(p); //p.getMass();

                    ModificationMatch modificationMatches[] = p.getVariableModifications();

                    ArrayList<String> mods = new ArrayList<>();
                    if(modificationMatches != null) {
                        for (ModificationMatch modificationMatch : modificationMatches) {
                            mods.add(modificationMatch.getSite() + ":" + modificationMatch.getModification());
                        }
                    }

                    // fixed modification
                    HashMap<Integer,Modification> fixedModificationMatches = ModificationGear.getInstance().getPossibleModificationSites(p,fixedPTMs);
                    for(Integer fix_i: fixedModificationMatches.keySet()){
                        mods.add(fix_i + ":" + fixedModificationMatches.get(fix_i).getName());
                    }


                    String modstring;
                    if(mods.size()>=1){
                        modstring = StringUtils.join(mods,";");
                    }else{
                        modstring = "-";
                    }

                    ByteArrayOutputStream bos = new ByteArrayOutputStream();

                    ObjectOutputStream oos = new ObjectOutputStream(bos);

                    oos.writeObject(p);

                    //pstmt.setString(1,p.getSequenceWithLowerCasePtms());
                    pstmt.setString(1,p.getSequence());
                    pstmt.setDouble(2,peptideMass);
                    pstmt.setBytes(3,bos.toByteArray());
                    pstmt.setString(4,proID);
                    pstmt.setInt(5,pos);
                    pstmt.setString(6,modstring);
                    pstmt.setString(7,p.getSequence());
                    pstmt.executeUpdate();

                }**/

            }


        }
        pstmt.close();

        reader.close();
        Cloger.getInstance().logger.info("Protein sequences:"+num);
        Cloger.getInstance().logger.info("Unique peptide sequences:"+searchedPeptides.size());

        connection.commit();
        stmt.executeUpdate("create index mass_index on prodb (mass)");
        stmt.executeUpdate("create index seq_index on prodb (sequence)");
        connection.commit();

        stmt.close();
        connection.close();


    }




}
