package main.java.pg;


import com.compomics.util.experiment.biology.enzymes.Enzyme;
import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.identification.protein_sequences.digestion.ExtendedPeptide;
import com.compomics.util.experiment.identification.protein_sequences.digestion.IteratorFactory;

import com.compomics.util.experiment.identification.protein_sequences.digestion.SequenceIterator;

import com.compomics.util.nucleotide.NucleotideSequenceImpl;

import com.compomics.util.parameters.identification.search.DigestionParameters;
import com.compomics.util.protein.AASequenceImpl;

import main.java.PSMMatch.HyperscoreMatch;
import main.java.util.Cloger;
import net.sf.jfasta.FASTAElement;
import net.sf.jfasta.FASTAFileReader;
import net.sf.jfasta.impl.FASTAElementIterator;
import net.sf.jfasta.impl.FASTAFileReaderImpl;
import org.apache.commons.cli.*;
import org.apache.commons.lang3.StringUtils;
import java.io.*;
import java.sql.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Taken DNA/RNA/VCF as input and pre-process the input to peptide sequences.
 */
public class InputProcessor {




    /**
     * 1 => protein
     * 2 => DNA
     * 3 => VCF
     * 4 => BED
     * 5 => GTF
     */
    public CParameter.TargetInputType inputType = CParameter.TargetInputType.protein;

    /**
     * The file
     */
    public String inputFile = "";

    public final static int InputTypeProtein = 1;
    public final static int InputTypeDNA = 2;
    public final static int InputTypeVCF = 3;
    public final static int InputTypeBED = 4;
    public final static int InputTypeGTF = 5;

    // The following two variables are used for VCF/BED/GTF input
    public String outdir = "./";
    public String annotationFolder = "./";

    // digest

    //public int minPepLength = 7;
    //public int maxPepLength = 45;

    /**
     * The minimum length of protein sequence to consider when translate DNA sequence to protein.
     * Default is 20.
     */
    public int minProLength = 20;

    /**
     * The frame to translate DNA sequence to protein
     * The right format is like this: "1,2,3,4,5,6","1,2,3","1".
     * "0" means to keep the longest frame. In default, for each frame only the longest protein is used.
     */
    public String frames = "0";

    /**
     * Protein reference database
     */
    public String referenceProDB = "";


    public static void main(String[] args) throws ParseException, IOException, SQLException, InterruptedException {


        Cloger.getInstance().logger.info("Start analysis");
        Cloger.getInstance().logger.info(StringUtils.join(args," "));

        HyperscoreMatch.generateFactorialValues(60);

        Options options = new Options();


        options.addOption("i", true, "Input file");
        options.addOption("db", true, "Fasta format database file");
        options.addOption("t",true,"Input type: 1=>protein,2=>DNA,3=>VCF,4=>BED,5=>GTF");
        options.addOption("f",true,"The frame to translate DNA sequence to protein. The right format is like this: \"1,2,3,4,5,6\",\"1,2,3\",\"1\". \"0\" means to keep the longest frame. In default, for each frame only the longest protein is used.");
        options.addOption("anno",true,"Annotation files folder for VCF/BED/GTF");

        options.addOption("h", false, "Help");


        CommandLineParser parser = new PosixParser();
        CommandLine cmd = parser.parse(options, args);
        if (cmd.hasOption("h") || cmd.hasOption("help") || args.length == 0) {
            HelpFormatter f = new HelpFormatter();
            f.printHelp("Options", options);
            System.exit(0);
        }
        String inputSeq = cmd.getOptionValue("i");
        String db = cmd.getOptionValue("db");
        InputProcessor inputProcessor = new InputProcessor();
        if(cmd.hasOption("t")){
            if(Integer.parseInt(cmd.getOptionValue("t"))==2){
                inputProcessor.setInputType(CParameter.TargetInputType.dna);
            }else if(Integer.parseInt(cmd.getOptionValue("t"))==3){
                inputProcessor.setInputType(CParameter.TargetInputType.vcf);
            }else if(Integer.parseInt(cmd.getOptionValue("t"))==4) {
                inputProcessor.setInputType(CParameter.TargetInputType.bed);
            }else if(Integer.parseInt(cmd.getOptionValue("t"))==5) {
                inputProcessor.setInputType(CParameter.TargetInputType.gtf);
            }
        }else{
            inputProcessor.setInputType(CParameter.TargetInputType.protein);
        }

        if(cmd.hasOption("f")){
            inputProcessor.frames = cmd.getOptionValue("f");
        }
        if(cmd.hasOption("anno")){
            inputProcessor.annotationFolder = cmd.getOptionValue("anno");
        }
        inputProcessor.referenceProDB = db;
        inputProcessor.run(inputSeq);
    }


    public  ArrayList<CPeptide> run(String inputSeq) throws SQLException, InterruptedException, IOException {
        ArrayList<CPeptide> cPeptides = new ArrayList<>();
        if(this.inputType.equals(CParameter.TargetInputType.dna)) {
            Cloger.getInstance().logger.info("The input is a DNA sequence.");
            cPeptides = translate(inputSeq);
        }else if(this.inputType.equals(CParameter.TargetInputType.protein)){
            // First we need to check whether the input is protein ID, or protein sequence or a
            // FASTA file including protein sequences

            File FI = new File(inputSeq);
            if(FI.isFile()) {
                Cloger.getInstance().logger.info("The input is a protein FASTA file.");
                cPeptides = digest(FI);
            }else{
                Cloger.getInstance().logger.info("The input is a protein sequence.");
                cPeptides = digest(inputSeq);
            }
        }else if(this.inputType.equals(CParameter.TargetInputType.vcf)){
            cPeptides = getSeqFromPGA(inputSeq);
        }else if(this.inputType.equals(CParameter.TargetInputType.bed)){
            cPeptides = getSeqFromPGA(inputSeq);
        }else if(this.inputType.equals(CParameter.TargetInputType.gtf)){
            cPeptides = getSeqFromPGA(inputSeq);
        }else{
            System.err.println("The input format is not supported:"+inputType);
            System.exit(0);
        }

        Cloger.getInstance().logger.info("Peptides:"+cPeptides.size());

        searchRefDB(cPeptides);

        for(CPeptide cPeptide: cPeptides){
            if(cPeptide.nhit==0) {
                // System.out.print("Peptides are not exists in the reference protein database: ");
                Cloger.getInstance().logger.info(cPeptide.peptideSequence+" is not exist in the reference protein database.");
                //System.out.println(cPeptide.peptideSequence + "\t" + cPeptide.nhit+"\t"+cPeptide.position+"\t"+cPeptide.frame+"\t"+cPeptide.proteinID);
            }
        }
        return(cPeptides);

    }

    private ArrayList<CPeptide> translate(String dnaSeq){

        ArrayList<CPeptide> cPeptides = new ArrayList<>();
        // Create a NucleotideSequenceImpl.
        NucleotideSequenceImpl nsi = new NucleotideSequenceImpl(dnaSeq);

        AASequenceImpl[] seqs = nsi.translate();
        if(this.frames.equalsIgnoreCase("0")) {
            Cloger.getInstance().logger.info("Frame for DNA sequence translation: longest ORF.");
            String longestProtein = "";
            int fm = -1;

            for (int i = 0; i < seqs.length; i++) {
                AASequenceImpl lSeq = seqs[i];
                String pep[] = lSeq.getSequence().split("_");
                // get the longest proteins.
                String tmp = Arrays.asList(pep).stream().max(Comparator.comparingInt(String::length)).get();
                if(tmp.length() > this.minProLength && tmp.length() > longestProtein.length()){
                    longestProtein = tmp;
                    fm = i+1;
                }
            }
            if(longestProtein.length() >= 1){
                try {
                    cPeptides = digest(longestProtein);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
                for(CPeptide cPeptide: cPeptides){
                    cPeptide.frame = fm;
                }
            }
        }else{
            Cloger.getInstance().logger.info("Frame for DNA sequence translation: "+this.frames+".");
            String fs[] = this.frames.split(",");
            HashSet<Integer> frs = new HashSet<>();
            for(int i=0;i<fs.length;i++){
                frs.add(Integer.valueOf(fs[i]));
            }
            for (int i = 0; i < seqs.length; i++) {
                int f = i+1;
                if(frs.contains(f)) {
                    AASequenceImpl lSeq = seqs[i];
                    String pep[] = lSeq.getSequence().split("_");
                    // get the longest proteins.
                    String tmp = Arrays.asList(pep).stream().max(Comparator.comparingInt(String::length)).get();
                    if (tmp.length() > this.minProLength) {
                        ArrayList<CPeptide> cps = null;
                        try {
                            cps = digest(tmp);
                        } catch (InterruptedException e) {
                            e.printStackTrace();
                        }
                        for(CPeptide cPeptide: cps){
                            cPeptide.frame = f;
                        }
                        cPeptides.addAll(cps);
                    }
                }
            }

        }

        return(cPeptides);
    }


    /**
     * Digest a protein into peptide sequences
     * @param proSeq A protein sequence or multiple protein sequences separated by ","
     * @return An HashMap<String,Integer> object which contains the peptides and their positions.
     */
    private ArrayList<CPeptide> digest(String proSeq) throws InterruptedException {
        ArrayList<CPeptide> cPeptides = new ArrayList<>();

        ArrayList<String> fixedModifications = new ArrayList<>();


        Enzyme enzyme = DatabaseInput.getEnzymeByIndex(CParameter.enzyme);
        DigestionParameters digestionPreferences = DatabaseInput.getDigestionPreferences(enzyme.getName(),CParameter.maxMissedCleavages);

        //DigestionPreferences digestionPreferences = DigestionPreferences.getDefaultPreferences();
        IteratorFactory iteratorModifications = new IteratorFactory(fixedModifications);

        ArrayList<String> proSeqs = new ArrayList<>();
        String protein_sequences[] = proSeq.split(",");
        Cloger.getInstance().logger.info("The number of protein sequences:"+protein_sequences.length);
        for(String protein_seq: protein_sequences) {

            if (protein_seq.contains("*")) {
                String seqs[] = protein_seq.split("\\*");
                for (int i = 0; i < seqs.length; i++) {
                    proSeqs.add(seqs[i]);
                }
            } else {
                proSeqs.add(protein_seq);
            }
        }

        HashMap<String,Integer> pepSet = new HashMap<>();

        for(String proteinSeq : proSeqs) {


            // digest
            SequenceIterator sequenceIterator = iteratorModifications.getSequenceIterator(proteinSeq, digestionPreferences, CParameter.minPeptideMass, CParameter.maxPeptideMass);

            ExtendedPeptide peptideWithPosition;
            while ((peptideWithPosition = sequenceIterator.getNextPeptide()) != null) {

                Peptide peptide = peptideWithPosition.peptide;
                if (peptide.getSequence().length() >= CParameter.minPeptideLength && peptide.getSequence().length() <= CParameter.maxPeptideLength) {
                    if(!pepSet.containsKey(peptide.getSequence())) {
                        // make sure a unique peptide sequence only has one copy in cPeptides
                        CPeptide cPeptide = new CPeptide();
                        cPeptide.peptideSequence = peptide.getSequence();
                        cPeptide.position = peptideWithPosition.position;
                        cPeptides.add(cPeptide);
                        pepSet.put(peptide.getSequence(),1);
                        //pep2pos.put(peptide.getSequence(),pos);
                        //System.out.println("\t"+peptide.getSequence()+"\t"+pos);
                    }
                }
            }
        }

        return(cPeptides);

    }


    /**
     * Digest novel protein sequences
     * @param dbFile Novel protein database file in fasta format
     * @return
     */
    private ArrayList<CPeptide> digest(File dbFile) throws IOException {

        HashMap<String,CPeptide> pepSet = new HashMap<>();
        //File dbFile = new File(fa);
        FASTAFileReader reader = new FASTAFileReaderImpl(dbFile);
        FASTAElementIterator it = reader.getIterator();
        int num = 0;
        while (it.hasNext()) {
            FASTAElement el = it.next();
            el.setLineLength(1);
            String headLine[] = el.getHeader().split("\\s+");
            String proID = headLine[0];
            num++;
            String proSeq = el.getSequence().toUpperCase();

            proSeq = proSeq.replaceAll("\\*$", "");

            ArrayList<CPeptide> cpeps = null;
            try {
                cpeps = digest(proSeq);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }


            for(CPeptide cPeptide:cpeps){
                cPeptide.proteinID.add(proID);

                if(pepSet.containsKey(cPeptide.peptideSequence)){
                    pepSet.get(cPeptide.peptideSequence).proteinID.add(proID);
                }else{
                    pepSet.put(cPeptide.peptideSequence,cPeptide);
                }
            }

            //proSeq = proSeq.replaceAll("I", "L");
        }
        reader.close();
        Cloger.getInstance().logger.info("Protein sequences:"+num);

        String out_pep_pro_file = this.outdir+"/pep2pro.txt";
        BufferedWriter pepWriter = new BufferedWriter(new FileWriter(new File(out_pep_pro_file)));
        Cloger.getInstance().logger.info("Output peptide to protein mapping information to file: "+out_pep_pro_file);
        pepWriter.write("peptide\tprotein\n");

        ArrayList<CPeptide> cPeptides = new ArrayList<>();
        for(String pep: pepSet.keySet()){
            cPeptides.add(pepSet.get(pep));
            pepWriter.write(pep+"\t"+StringUtils.join(pepSet.get(pep).proteinID,";")+"\n");
        }
        pepWriter.close();
        Cloger.getInstance().logger.info("Peptide sequences:"+pepSet.size());
        return(cPeptides);

    }


    public void searchRefDB(ArrayList<CPeptide> cPeptides) throws SQLException {
        if(!this.referenceProDB.isEmpty()){
            String db = referenceProDB;
            // Search spectra matched to a target peptide against a reference protein database.
            // Consider both fixed and variable modifications.
            if(referenceProDB.toLowerCase().endsWith(".fa") || referenceProDB.toLowerCase().endsWith(".fasta")) {

                // For indexed database ended with ".sqldb"
                String sqldb_file = referenceProDB+".sqldb";
                File SQLDB = new File(sqldb_file);
                // use sql database
                if(SQLDB.isFile()){
                    db = sqldb_file;
                    Cloger.getInstance().logger.info("Use indexed database:"+sqldb_file);

                    // read digested peptides from SQL database
                    Connection connection = DriverManager.getConnection("jdbc:sqlite:" + db);
                    PreparedStatement pstmt = connection.prepareStatement("select 1 from prodb where sequence = ? limit 1;");
                    for(CPeptide cPeptide : cPeptides){
                        //System.out.println(cPeptide.peptideSequence);
                        pstmt.setString(1,cPeptide.peptideSequence.replaceAll("I","L"));
                        ResultSet res = pstmt.executeQuery();
                        int count = 0 ;
                        while (res.next()){
                            count = res.getInt(1);
                        }
                        cPeptide.nhit = count;
                    }
                    pstmt.close();
                    connection.close();
                    //System.out.println("Peptides from reference database:"+npep);

                }else {
                    // No indexed database
                    Cloger.getInstance().logger.info("Don't find indexed database:"+sqldb_file);
                    Cloger.getInstance().logger.info("Use database:"+referenceProDB);

                    if(CParameter.peptideMappingType == CParameter.PepMappingType_FM) {

                        PepMapping pepMapping = PepMapping.getInstance();
                        pepMapping.loadDB(referenceProDB);
                        for (CPeptide cPeptide : cPeptides) {
                            if(pepMapping.hasProtein(cPeptide.peptideSequence)){
                                cPeptide.nhit = 1;
                            } else {
                                cPeptide.nhit = 0;
                            }

                        }
                        //pepMapping.close();
                    }else{
                        // digest
                        HashSet<String> pepSet = buildPeptideDB(db);

                        for (CPeptide cPeptide : cPeptides) {
                            if (pepSet.contains(cPeptide.peptideSequence.replaceAll("I", "L"))) {
                                cPeptide.nhit = 1;
                            } else {
                                cPeptide.nhit = 0;
                            }

                        }
                    }
                }
            }else{
                System.err.println("Please provide valid protein database: .fa or .fasta format.");
                System.exit(0);
            }


        }
    }


    /**
     * Digest proteins
     * @param db protein database
     * @return HashSet<String> peptide hash set
     */
    private HashSet<String> buildPeptideDB(String db){

        long startTime=System.currentTimeMillis();

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // read protein database
        File dbFile = new File(db);
        FASTAFileReader reader = null;
        try {
            reader = new FASTAFileReaderImpl(dbFile);
        } catch (IOException e) {
            e.printStackTrace();
        }
        FASTAElementIterator it = null;
        try {
            it = reader.getIterator();
        } catch (IOException e) {
            e.printStackTrace();
        }

        // digest protein
        Enzyme enzyme = DatabaseInput.getEnzymeByIndex(CParameter.enzyme);
        DigestionParameters digestionPreferences = DatabaseInput.getDigestionPreferences(enzyme.getName(),CParameter.maxMissedCleavages);

        Cloger.getInstance().logger.info("Enzyme:"+enzyme.getName()+", max missed cleavages:"+ CParameter.maxMissedCleavages);

        //HashSet<String> searchedPeptides = new HashSet<>();


        int cpu = Runtime.getRuntime().availableProcessors();
        ExecutorService fixedThreadPool = Executors.newFixedThreadPool(cpu);

        ConcurrentHashMap<String, HashSet<String>> pro2pep = new ConcurrentHashMap<>();

        int num = 0;
        try {
            while (it.hasNext()) {
                FASTAElement el = it.next();
                el.setLineLength(1);
                String headLine[] = el.getHeader().split("\\s+");
                String proID = headLine[0];
                num++;
                String proSeq = el.getSequence();
                // digest
                fixedThreadPool.execute(new PDigestProteinWorker(proID, proSeq, digestionPreferences, pro2pep));

            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        fixedThreadPool.shutdown();

        try {
            fixedThreadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        reader.close();
        Cloger.getInstance().logger.info("Protein sequences:" + num);


        // peptides from reference database
        HashSet<String> pepSet = new HashSet<>();

        for (String proID : pro2pep.keySet()) {
            pepSet.addAll(pro2pep.get(proID));
        }

        long bTime = System.currentTimeMillis();
        Cloger.getInstance().logger.info("Total unique peptides: "+pepSet.size());
        Cloger.getInstance().logger.info("Time used for building peptide index:" + (bTime - startTime) / 1000 + " s.");


        return(pepSet);
    }



    public void setInputType(CParameter.TargetInputType inputType) {
        this.inputType = inputType;
    }


    /**
     * Generate peptide sequences based on VCF/BED/GTF file.
     * @param input VCF/BED/GTF file
     * @throws InterruptedException
     */
    public ArrayList<CPeptide> getSeqFromPGA(String input) throws InterruptedException, IOException {
        String rcode = "";
        String outfa = "";
        ArrayList<CPeptide> cPeptides = new ArrayList<>();

        // For VCF/BED/GTF, there are only two options.
        String bool_get_longest = "FALSE";
        if(this.frames.equalsIgnoreCase("0")){
            // "0" means to keep the longest frame. In default, for each frame only the longest protein is used.
            bool_get_longest = "TRUE";
        }

        // the file extension is very important
        if(input.toLowerCase().endsWith(".vcf")){
            rcode = "library(\"PGA\");library(BSgenome.Hsapiens.UCSC.hg19);\n" +
                    "dbfile <- dbCreator(vcfFile=\"" + input + "\"," +
                    "                    annotation_path=\"" + annotationFolder + "\"," +
                    "                    outdir=\"" + outdir + "\"," +
                    "                    make_decoy=FALSE," +
                    "                    bool_get_longest=" + bool_get_longest +","+
                    "                    outfile_name=\"pga\"" + "," +
                    "                    genome=Hsapiens);";

            runR(rcode);

            outfa = this.outdir + "/pga_snv.fasta";
        }else if(input.toLowerCase().endsWith(".bed")){
            rcode = "library(\"PGA\");library(BSgenome.Hsapiens.UCSC.hg19);\n" +
                    "dbfile <- dbCreator(bedFile=\"" + input + "\",\n" +
                    "                    annotation_path=\"" + annotationFolder + "\"," +
                    "                    outdir=\"" + outdir + "\"," +
                    "                    make_decoy=FALSE," +
                    "                    bool_get_longest=" + bool_get_longest +","+
                    "                    outfile_name=\"pga\"" + "," +
                    "                    genome=Hsapiens);";
            runR(rcode);

            outfa = this.outdir + "/pga_junc.fasta";
        }else if(input.toLowerCase().endsWith(".gtf")){
            rcode = "library(\"PGA\");library(BSgenome.Hsapiens.UCSC.hg19);\n" +
                    "dbfile <- dbCreator(gtfFile=\"" + input + "\",\n" +
                    "                    annotation_path=\"" + annotationFolder + "\"," +
                    "                    outdir=\"" + outdir + "\"," +
                    "                    make_decoy=FALSE," +
                    "                    bool_get_longest=" + bool_get_longest +","+
                    "                    outfile_name=\"pga\"" + "," +
                    "                    genome=Hsapiens);";
            runR(rcode);

            outfa = this.outdir + "/pga_ntx.fasta";
        }else{
            System.err.println("The file format is not supported:"+input);
            System.exit(0);
        }

        // read the FASTA file
        File dbFile = new File(outfa);
        FASTAFileReader reader = new FASTAFileReaderImpl(dbFile);
        FASTAElementIterator it = reader.getIterator();

        // remove redundancy
        HashMap<String,CPeptide> cPeptideHashMap = new HashMap<>();

        while (it.hasNext()) {
            FASTAElement el = it.next();
            el.setLineLength(1);
            String headLine[] = el.getHeader().split("\\s+");
            String proID = headLine[0];
            String proSeq = el.getSequence().toUpperCase();

            // SNV fasta header
            // >SNV11|ENSP00000360695|S86L
            proSeq = proSeq.replaceAll("\\*","");

            // Currently, if the input is VCF, then we need to use the peptides that covers the SNV site.
            int snvPos = -1;
            if(this.inputType.equals(CParameter.TargetInputType.vcf)) {
                String[] ids = proID.split("\\|");
                snvPos = Integer.valueOf(ids[ids.length - 1].replaceAll("[^0-9]", ""));
                Cloger.getInstance().logger.info(proID + "\t" + snvPos);
            }

            ArrayList<CPeptide> cPeps = digest(proSeq);
            //cPeps.stream().forEach(cPeptide -> cPeptide.proteinID = proID);
            for(CPeptide cp : cPeps){
                cp.proteinID.add(proID);
            }
            for(CPeptide cPeptide : cPeps){
                if(this.inputType.equals(CParameter.TargetInputType.vcf)){
                    if(snvPos < cPeptide.position || snvPos > (cPeptide.position+cPeptide.peptideSequence.length()) ) {
                        continue;
                    }
                }
                if(cPeptideHashMap.containsKey(cPeptide.peptideSequence)){
                    cPeptideHashMap.get(cPeptide.peptideSequence).proteinID.add(proID);
                    //cPeptideHashMap.get(cPeptide.peptideSequence).proteinID = cPeptideHashMap.get(cPeptide.peptideSequence).proteinID + ";" + proID;
                }else{
                    cPeptideHashMap.put(cPeptide.peptideSequence,cPeptide);
                }
            }

        }
        reader.close();
        cPeptides.addAll(cPeptideHashMap.values());
        return(cPeptides);
    }


    /**
     * Run R code
     * @param code R code
     * @throws InterruptedException
     */
    public void runR(String code) throws InterruptedException {

        try {

            String rbin= "R";

            Process p = null;
            if(System.getProperty("os.name").equals("Linux")){
                rbin= ""+rbin+"";
                rbin = rbin + " --vanilla --slave";
                p = Runtime.getRuntime().exec(rbin);
            }else if (System.getProperty("os.name").startsWith("Windows")) {
                rbin= "\""+rbin+"\"";
                rbin = rbin + " --vanilla --slave";
                rbin = "cmd /c "+rbin;
                Cloger.getInstance().logger.info(rbin);
                p = Runtime.getRuntime().exec(rbin);
            }else{
                rbin= ""+rbin+"";
                rbin = rbin + " --vanilla --slave";
                p = Runtime.getRuntime().exec(rbin);
            }


            BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
            BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(p.getOutputStream()));

            BufferedInputStream errin = new BufferedInputStream(p.getErrorStream());
            BufferedReader errbr = new BufferedReader(new InputStreamReader(errin));
            Cloger.getInstance().logger.info(code);
            writer.write(code + "\n");
            writer.flush();
            writer.close();
            String line;
            p.waitFor();
            if (true) {
                while ((line = reader.readLine()) != null) {
                    Cloger.getInstance().logger.info(line); // reading from process
                }
            }
            if (true) {
                while ((line = errbr.readLine()) != null) {
                    Cloger.getInstance().logger.info(line); // reading from process
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

}
