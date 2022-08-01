package main.java.pg;

import com.compomics.util.experiment.biology.enzymes.Enzyme;
import com.compomics.util.experiment.biology.enzymes.EnzymeFactory;
import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.identification.protein_sequences.digestion.ExtendedPeptide;
import com.compomics.util.experiment.identification.protein_sequences.digestion.IteratorFactory;

import com.compomics.util.experiment.identification.protein_sequences.digestion.SequenceIterator;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;

import com.compomics.util.parameters.identification.search.DigestionParameters;

import com.compomics.util.pride.CvTerm;
import main.java.OpenModificationSearch.DigestProteinWorker;
import main.java.OpenModificationSearch.GeneratePeptideformWorker;
import main.java.OpenModificationSearch.ModificationDB;
import main.java.PSMMatch.JPeptide;
import main.java.PSMMatch.ScoreResult;
import main.java.plot.PeakAnnotation;
import main.java.util.Cloger;
import net.sf.jfasta.FASTAElement;
import net.sf.jfasta.FASTAFileReader;
import net.sf.jfasta.impl.FASTAElementIterator;
import net.sf.jfasta.impl.FASTAFileReaderImpl;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;

import java.io.*;
import java.sql.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * search an MS/MS spectrum against a reference protein database
 */
public final class DatabaseInput {

    public String db = null;

    public HashSet<String> target_peptides = new HashSet<>();

    public void add_target_peptides(ArrayList<String> peps){
        for(String pep : peps) {
            pep.replaceAll("I","L");
            this.target_peptides.add(pep.replaceAll("I","L"));
        }
    }

    /**
     * fixed modification list: ptm name
     */
    public ArrayList<String> fixedModifications = new ArrayList<>();


    /**
     * variable modification list: ptm name
     */
    public ArrayList<String> varModifications = new ArrayList<String>();


    public void addFixedModifications(int id) {
        if (!fixedModifications.contains(ModificationGear.getInstance().getPTMname(id))) {
            fixedModifications.add(ModificationGear.getInstance().getPTMname(id));
        }
    }

    /**
     * A file contains tags for spectra
     */
    public String tagFile = "";


    /**
     * add variable modification according to the int id of ptm
     * @param id
     */
    public void addVarModifications(int id) {
        if (!varModifications.contains(ModificationGear.getInstance().getPTMname(id))) {
            varModifications.add(ModificationGear.getInstance().getPTMname(id));
        }
    }

    /**
     * saved the spectra that will need to map to peptides
     */
    public HashMap<String, Spectrum> mSnSpectrums = new HashMap<>();

    /**
     * key: spectrum title, value: an ArrayList of Peptide objects
     */
    public ConcurrentHashMap<String,ArrayList<JPeptide>> ms2peptide = new ConcurrentHashMap<>();



    public DatabaseInput(String file){
        this.db = file;

    }

    /**
     * Set the enzyme and missed cleavage
     * @param enzymeName Enzyme name
     * @param enzymeMissedCleavages The max missed cleavage
     * @return
     */
    public static DigestionParameters getDigestionPreferences(String enzymeName, int enzymeMissedCleavages) {
        DigestionParameters digestionParameters = new DigestionParameters();
        if(EnzymeFactory.getInstance().enzymeLoaded(enzymeName)){
            digestionParameters.setCleavageParameter(DigestionParameters.CleavageParameter.enzyme);
            Enzyme enzyme = EnzymeFactory.getInstance().getEnzyme(enzymeName);
            digestionParameters.addEnzyme(enzyme);
            digestionParameters.setnMissedCleavages(enzymeName, enzymeMissedCleavages);
        }else if(enzymeName.equalsIgnoreCase("NoEnzyme")){
            digestionParameters.setCleavageParameter(DigestionParameters.CleavageParameter.unSpecific);
            Enzyme enzyme = getEnzymeByIndex(0);
            digestionParameters.addEnzyme(enzyme);
            digestionParameters.setnMissedCleavages(enzymeName, 1000);
        }else{
            System.err.println("No valid enzyme:"+enzymeName);
            System.exit(0);
        }

        return digestionParameters;
    }

    /**
     * In default, using trypsin as the enzyme with up to 2 missed cleavage.
     * @return
     */
    public static DigestionParameters getDigestionPreferences(){
        return(getDigestionPreferences("Trypsin",2));
    }


    /**
     * Get the enzyme object according to enzyme index
     * @param ind enzyme index, 0-based index.
     * @return An enzyme object
     */
    public static Enzyme getEnzymeByIndex(int ind){

        ArrayList<Enzyme> enzymes = new ArrayList<>();

        // 0 non-specific digestion
        Enzyme enzyme = new Enzyme("NoEnzyme");
        // no useful
        enzyme.addAminoAcidBefore('R');
        // no useful
        enzyme.addAminoAcidBefore('K');
        enzyme.setCvTerm(new CvTerm("PSI-MS", "MS:1001956", "NoEnzyme", null));
        enzymes.add(enzyme);

        enzyme = new Enzyme("Trypsin");
        enzyme.addAminoAcidBefore('R');
        enzyme.addAminoAcidBefore('K');
        enzyme.addRestrictionAfter('P');
        enzyme.setCvTerm(new CvTerm("PSI-MS", "MS:1001251", "Trypsin", null));
        enzymes.add(enzyme);

        enzyme = new Enzyme("Trypsin (no P rule)");
        enzyme.addAminoAcidBefore('R');
        enzyme.addAminoAcidBefore('K');
        enzyme.setCvTerm(new CvTerm("PSI-MS", "MS:1001313", "Trypsin/P", null));
        enzymes.add(enzyme);

        enzyme = new Enzyme("Arg-C");
        enzyme.addAminoAcidBefore('R');
        enzyme.addRestrictionAfter('P');
        enzyme.setCvTerm(new CvTerm("PSI-MS", "MS:1001303", "Arg-C", null));
        enzymes.add(enzyme);

        enzyme = new Enzyme("Arg-C (no P rule)");
        enzyme.addAminoAcidBefore('R');
        enzymes.add(enzyme);

        enzyme = new Enzyme("Arg-N");
        enzyme.addAminoAcidAfter('R');
        enzymes.add(enzyme);

        enzyme = new Enzyme("Glu-C");
        enzyme.addAminoAcidBefore('E');
        enzyme.setCvTerm(new CvTerm("PSI-MS", "MS:1001917", "glutamyl endopeptidase", null));
        enzymes.add(enzyme);

        enzyme = new Enzyme("Lys-C");
        enzyme.addAminoAcidBefore('K');
        enzyme.addRestrictionAfter('P');
        enzyme.setCvTerm(new CvTerm("PSI-MS", "MS:1001309", "Lys-C", null));
        enzymes.add(enzyme);

        enzyme = new Enzyme("Lys-C (no P rule)");
        enzyme.addAminoAcidBefore('K');
        enzyme.setCvTerm(new CvTerm("PSI-MS", "MS:1001310", "Lys-C/P", null));
        enzymes.add(enzyme);

        enzyme = new Enzyme("Lys-N");
        enzyme.addAminoAcidAfter('K');
        enzymes.add(enzyme);

        enzyme = new Enzyme("Asp-N");
        enzyme.addAminoAcidAfter('D');
        enzyme.setCvTerm(new CvTerm("PSI-MS", "MS:1001304", "Asp-N", null));
        enzymes.add(enzyme);

        enzyme = new Enzyme("Asp-N (ambic)");
        enzyme.addAminoAcidAfter('D');
        enzyme.addAminoAcidAfter('E');
        enzyme.setCvTerm(new CvTerm("PSI-MS", "MS:1001305", "Asp-N_ambic", null));
        enzymes.add(enzyme);

        enzyme = new Enzyme("Chymotrypsin");
        enzyme.addAminoAcidBefore('F');
        enzyme.addAminoAcidBefore('Y');
        enzyme.addAminoAcidBefore('W');
        enzyme.addAminoAcidBefore('L');
        enzyme.addRestrictionAfter('P');
        enzyme.setCvTerm(new CvTerm("PSI-MS", "MS:1001306", "Chymotrypsin", null));
        enzymes.add(enzyme);

        enzyme = new Enzyme("Chymotrypsin (no P rule)");
        enzyme.addAminoAcidBefore('F');
        enzyme.addAminoAcidBefore('Y');
        enzyme.addAminoAcidBefore('W');
        enzyme.addAminoAcidBefore('L');
        enzymes.add(enzyme);

        enzyme = new Enzyme("Pepsin A");
        enzyme.addAminoAcidBefore('F');
        enzyme.addAminoAcidBefore('L');
        enzyme.setCvTerm(new CvTerm("PSI-MS", "MS:1001311", "Pepsin A", null));
        enzymes.add(enzyme);

        enzyme = new Enzyme("CNBr");
        enzyme.addAminoAcidBefore('M');
        enzyme.setCvTerm(new CvTerm("PSI-MS", "MS:1001307", "CNBr", null));
        enzymes.add(enzyme);

        enzyme = new Enzyme("Thermolysin");
        enzyme.addAminoAcidAfter('A');
        enzyme.addAminoAcidAfter('F');
        enzyme.addAminoAcidAfter('I');
        enzyme.addAminoAcidAfter('L');
        enzyme.addAminoAcidAfter('M');
        enzyme.addAminoAcidAfter('V');
        enzymes.add(enzyme);

        enzyme = new Enzyme("LysargiNase");
        enzyme.addAminoAcidAfter('R');
        enzyme.addAminoAcidAfter('K');
        enzymes.add(enzyme);

        //if(ind<=0 || ind > enzymes.size()){
        if(ind < 0 || ind > enzymes.size()){
            System.err.println("Please provide a valid enzyme number:"+ind);
            System.exit(0);
        }
        // System.out.println("Use enzyme:"+enzymes.get(ind-1).getName());
        Cloger.getInstance().logger.info("Use enzyme:"+enzymes.get(ind).getName());

        // return(enzymes.get(ind-1));
        return(enzymes.get(ind));
    }


    /**
     * Extract reference peptide matches for each target spectrum.
     * @throws IOException
     * @throws InterruptedException
     */
    public void readDB() throws IOException, InterruptedException {
        File dbFile = new File(this.db);
        FASTAFileReader reader = new FASTAFileReaderImpl(dbFile);
        FASTAElementIterator it = reader.getIterator();

        // variable modifications
        ModificationFactory ptmFactory = ModificationFactory.getInstance();
        /**
        ArrayList<Modification> varPTMs = new ArrayList<>();
        for(String ptmName: varModifications){
            varPTMs.add(ptmFactory.getModification(ptmName));
        }**/

        // digest protein
        //DigestionPreferences digestionPreferences = DigestionPreferences.getDefaultPreferences();
        Enzyme enzyme = getEnzymeByIndex(CParameter.enzyme);
        DigestionParameters digestionPreferences = getDigestionPreferences(enzyme.getName(),CParameter.maxMissedCleavages);

        ArrayList<String> fixM = new ArrayList<>();
        IteratorFactory iteratorModifications = new IteratorFactory(fixM);

        Cloger.getInstance().logger.info("Enzyme: "+enzyme.getName()+", maxMissedCleavages: "+CParameter.maxMissedCleavages);

        // prepare spectrum index for fast peptide query
        HashMap<Integer,ArrayList<String>> spectraIndex = new HashMap<>();
        // step: 0.1
        for(String title : SpectraInput.spectraMap.keySet()){
            Spectrum mSnSpectrum = SpectraInput.spectraMap.get(title);
            double mass = mSnSpectrum.getPrecursor().getMass(mSnSpectrum.getPrecursor().possibleCharges[0]);
            int va = (int) Math.round(10.0*mass);
            if(spectraIndex.containsKey(va)){
                spectraIndex.get(va).add(mSnSpectrum.getSpectrumTitle());
            }else{
                ArrayList<String> spectra_titles = new ArrayList<>();
                spectra_titles.add(mSnSpectrum.getSpectrumTitle());
                spectraIndex.put(va,spectra_titles);
            }
        }

        // if a peptide is searched before, then it will be omitted.
        ConcurrentHashMap<String,Integer> searchedPeptides = new ConcurrentHashMap<>();

        ExecutorService fixedThreadPool = Executors.newFixedThreadPool(CParameter.cpu);
        /**
         * key: spectrum title, value: an ArrayList of Peptide objects
         */
        ConcurrentHashMap<String,HashMap<String,ArrayList<Peptide>>> protein2ms2peptide = new ConcurrentHashMap<>();

        int total_proteins = 0;
        while (it.hasNext()) {
            FASTAElement el = it.next();
            el.setLineLength(1);
            String headLine[] = el.getHeader().split("\\s+");
            String proID = headLine[0];
            total_proteins++;
            // System.out.println(proID);
            String proSeq = el.getSequence().toUpperCase();
            protein2ms2peptide.put(proID,new HashMap<>());

            if(proSeq.length() < CParameter.minPeptideLength){
                continue;
            }

            proSeq = proSeq.replaceAll("\\*","");

            fixedThreadPool.execute(new ProteinToSpectrumMatch(iteratorModifications,
                    proSeq,
                    digestionPreferences,
                    searchedPeptides,
                    spectraIndex,
                    protein2ms2peptide.get(proID)));

            //fixedThreadPool.execute(new DBdigestProteinWorker(spectraIndex, proID, proSeq, digestionPreferences, ms2peptide, digestedPeptides));
            //SequenceIterator sequenceIterator = iteratorModifications.getSequenceIterator(proSeq, digestionPreferences, 600.0, 5000.0);

        }
        reader.close();

        fixedThreadPool.shutdown();
        fixedThreadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);

        HashMap<String, HashSet<String>> spectrum2peptideform = new HashMap<>();

        for(String proID: protein2ms2peptide.keySet()){
            for(String spectrum_id: protein2ms2peptide.get(proID).keySet()){
                for(Peptide peptide: protein2ms2peptide.get(proID).get(spectrum_id)){

                    if(CParameter.search_type.equals(CParameter.SearchType.known)){
                        if(this.target_peptides.contains(peptide.getSequence())){
                            continue;
                        }
                    }

                    String pep_mod = peptide.getSequence()+"|"+ModificationDB.getInstance().getModificationString(peptide);
                    if(spectrum2peptideform.containsKey(spectrum_id)){
                        if(spectrum2peptideform.get(spectrum_id).contains(pep_mod)){
                            continue;
                        }else{
                            if(ms2peptide.containsKey(spectrum_id)) {
                                ms2peptide.get(spectrum_id).add(new JPeptide(peptide));
                            }else{
                                ms2peptide.put(spectrum_id,new ArrayList<>());
                                ms2peptide.get(spectrum_id).add(new JPeptide(peptide));
                            }
                            spectrum2peptideform.get(spectrum_id).add(pep_mod);
                        }
                    }else{
                        if(ms2peptide.containsKey(spectrum_id)) {
                            ms2peptide.get(spectrum_id).add(new JPeptide(peptide));
                        }else{
                            ms2peptide.put(spectrum_id,new ArrayList<>());
                            ms2peptide.get(spectrum_id).add(new JPeptide(peptide));
                        }
                        spectrum2peptideform.put(spectrum_id,new HashSet<>());
                        spectrum2peptideform.get(spectrum_id).add(pep_mod);
                    }
                }
            }
        }

        spectraIndex = null;
        Cloger.getInstance().logger.info("Protein sequences:"+total_proteins);

    }



    public static ConcurrentHashMap<String,ArrayList<JPeptide>> protein_digest(String db){

        long startTime=System.currentTimeMillis();
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
        DigestionParameters digestionParameters = DatabaseInput.getDigestionPreferences(enzyme.getName(),CParameter.maxMissedCleavages);

        // protein digestion
        int cpu = 1;
        if(CParameter.cpu <= 0){
            cpu = Runtime.getRuntime().availableProcessors();
        }else{
            cpu = CParameter.cpu;
        }
        ExecutorService fixedThreadPool = Executors.newFixedThreadPool(cpu);

        // ConcurrentHashMap<String, HashMap<String, Double>> pro2pep = new ConcurrentHashMap<>();
        // ConcurrentHashMap<String, ArrayList<Peptide>> pro2pep = new ConcurrentHashMap<>();
        ConcurrentHashMap<String, HashSet<String>> pro2pepSeq = new ConcurrentHashMap<>(20000,0.75F, cpu);

        int num = 0;
        try {
            while (it.hasNext()) {
                FASTAElement el = it.next();
                el.setLineLength(1);
                String headLine[] = el.getHeader().split("\\s+");
                String proID = headLine[0];
                num++;
                String proSeq = el.getSequence();
                fixedThreadPool.execute(new DigestProteinWorker(proID, proSeq, digestionParameters, pro2pepSeq, true));
                //fixedThreadPool.execute(new DigestProteinWorker(proID, proSeq, digestionParameters, pro2pep));

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

        ConcurrentHashMap<String,ArrayList<JPeptide>> searchedPeptides = new ConcurrentHashMap<>(10000, 0.75F, cpu);
        for (String proID : pro2pepSeq.keySet()) {
            for (String pep : pro2pepSeq.get(proID)){
                searchedPeptides.put(pep, new ArrayList<>());
            }
        }
        pro2pepSeq.clear();
        pro2pepSeq = null;

        long bTime = System.currentTimeMillis();
        Cloger.getInstance().logger.info("Protein sequences:" + num+", total unique peptide sequences:"+searchedPeptides.size());
        Cloger.getInstance().logger.info("Time used for protein digestion:" + (bTime - startTime) / 1000 + " s.");
        // protein digestion done.
        return searchedPeptides;
    }

    public static ConcurrentHashMap<String,ArrayList<JPeptide>> protein_digest(String db, HashSet<String> target_peptides){
        ConcurrentHashMap<String,ArrayList<JPeptide>> searchedPeptides = protein_digest(db);
        for(String pep: target_peptides){
            if(searchedPeptides.containsKey(pep)){
                searchedPeptides.remove(pep);
            }
        }
        return searchedPeptides;
    }

    public static HashMap<Long, ArrayList<JPeptide>> generate_reference_peptide_index(String db, HashSet<String> target_peptides, double mass_resolution){

        long startTime=System.currentTimeMillis();
        ConcurrentHashMap<String,ArrayList<JPeptide>> searchedPeptides = protein_digest(db,target_peptides);

        GeneratePeptideformWorker.set_mass_range(SpectraInput.spectraMap);

        int cpu = 1;
        if(CParameter.cpu <= 0){
            cpu = Runtime.getRuntime().availableProcessors();
        }else{
            cpu = CParameter.cpu;
        }
        ExecutorService fixedThreadPool = Executors.newFixedThreadPool(cpu);
        // generate peptide forms for all unique peptide sequences.
        for(String peptide_seq: searchedPeptides.keySet()){
            if(!target_peptides.contains(peptide_seq)) {
                fixedThreadPool.execute(new GeneratePeptideformWorker(peptide_seq, searchedPeptides));
            }
        }

        fixedThreadPool.shutdown();

        try {
            fixedThreadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        GeneratePeptideformWorker.reset_mass_range();

        // generate peptide forms for all unique peptide sequences finished.
        Cloger.getInstance().logger.info("Generating peptide forms is done!");

        // generate peptide indexing for fast unrestricted modified peptide matching
        long n_peptide_forms = 0;
        HashMap<Long,ArrayList<JPeptide>> peptideIndexMap = new HashMap<>(10000);
        for (String peptide_seq : searchedPeptides.keySet()) {
            for (JPeptide pep : searchedPeptides.get(peptide_seq)) { //at main.java.OpenModificationSearch.ModificationDB.buildPeptideDB(ModificationDB.java:874) Exception in thread "main" java.lang.NullPointerException
                n_peptide_forms = n_peptide_forms + 1;
                long intValue = Math.round(pep.getMass()*mass_resolution);
                if (!peptideIndexMap.containsKey(intValue)) {
                    peptideIndexMap.put(intValue, new ArrayList<>(1000));
                }
                peptideIndexMap.get(intValue).add(pep);
            }
        }

        // sort peptide based on peptide mass: from small to large
        Cloger.getInstance().logger.info("Sort peptide index ...");
        for(long intValue: peptideIndexMap.keySet()){
            Collections.sort(peptideIndexMap.get(intValue), comparator_peptide_mass);
        }
        Cloger.getInstance().logger.info("Sort peptide index done");


        long bTime = System.currentTimeMillis();

        Cloger.getInstance().logger.info("Total unique peptide: "+searchedPeptides.size() + ", total peptide forms:"+n_peptide_forms);


        searchedPeptides.clear();
        searchedPeptides = null;
        Cloger.getInstance().logger.info("Index block size:"+peptideIndexMap.size());
        Cloger.getInstance().logger.info("Time used for building peptide index:" + (bTime - startTime) / 1000 + " s.");
        return peptideIndexMap;
    }

    /**
     * Sort an array of JPeptide objects based on peptide mass: from small to large.
     */
    private static final Comparator<JPeptide> comparator_peptide_mass = (s1, s2) -> {
        if (s1.getMass() < s2.getMass()) {
            return -1;
        } else if (s1.getMass() > s2.getMass()) {
            return 1;
        } else {
            return 0;
        }
    };


    public HashMap<String, Double> get_spectra2score(String psm_file, String score_direction){
        // read PSM file and get the lowest score for each spectrum
        BufferedReader bReader = null;
        try {
            bReader = new BufferedReader(new FileReader(new File(psm_file)));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        String line = null;
        try {
            line = bReader.readLine().trim();
        } catch (IOException e) {
            e.printStackTrace();
        }
        String []headLine = line.split("\t");
        HashMap<String,Integer> headMap = new HashMap<>();
        for(int i=0;i<headLine.length;i++){
            headMap.put(headLine[i],i);
        }

        // store the lowest score (min score) for each spectrum.
        HashMap<String, Double> spectra2score = new HashMap<>();

        try {
            while((line = bReader.readLine())!=null){
                line = line.trim();
                String []d = line.split("\t");
                String spectrum_title  = d[headMap.get("spectrum_title")];
                //String peptideSequence = d[headMap.get("peptide")];
                double score = Double.parseDouble(d[headMap.get("score")]);


                if(score_direction.equalsIgnoreCase("min")) {
                    // save the min score for each spectrum needed to validate.
                    if (spectra2score.containsKey(spectrum_title)) {
                        if (spectra2score.get(spectrum_title) > score) {
                            spectra2score.put(spectrum_title, score);
                        }
                    } else {
                        spectra2score.put(spectrum_title, score);
                    }
                }else{
                    // save the max score for each spectrum needed to validate.
                    if (spectra2score.containsKey(spectrum_title)) {
                        if (spectra2score.get(spectrum_title) < score) {
                            spectra2score.put(spectrum_title, score);
                        }
                    } else {
                        spectra2score.put(spectrum_title, score);
                    }
                }

            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        try {
            bReader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return spectra2score;
    }

    public HashMap<String, Double> get_spectra2score(ArrayList<PeptideInput> peptideInputs, String score_direction){
        // store the lowest score (min score) or max score for each spectrum.
        HashMap<String, Double> spectra2score = new HashMap<>();
        for (PeptideInput peptideInput : peptideInputs) {
            for (JPeptide jPeptide : peptideInput.getPtmIsoforms()) {
                for (int i = 0; i < jPeptide.spectraIndexs.size(); i++) {
                    String spectrum_title = jPeptide.spectraIndexs.get(i);
                    ScoreResult scoreResult = jPeptide.scores.get(i);
                    double score = scoreResult.score;
                    if(score_direction.equalsIgnoreCase("min")) {
                        // save the min score for each spectrum needed to validate.
                        if (spectra2score.containsKey(spectrum_title)) {
                            if (spectra2score.get(spectrum_title) > score) {
                                spectra2score.put(spectrum_title, score);
                            }
                        } else {
                            spectra2score.put(spectrum_title, score);
                        }
                    }else{
                        // save the max score for each spectrum needed to validate.
                        if (spectra2score.containsKey(spectrum_title)) {
                            if (spectra2score.get(spectrum_title) < score) {
                                spectra2score.put(spectrum_title, score);
                            }
                        } else {
                            spectra2score.put(spectrum_title, score);
                        }
                    }
                }
            }
        }
        return spectra2score;
    }


    public void db_search(ArrayList<PeptideInput> peptideInputs) throws IOException {

        HashMap<String, Double> spectra2score = get_spectra2score(peptideInputs,"max");

        long startTime=System.currentTimeMillis();
        HashMap<Long, ArrayList<JPeptide>> peptideIndexMap = generate_reference_peptide_index(this.db,this.target_peptides, 10.0);
        ModificationDB.getInstance().add_peptide_index(peptideIndexMap);

        int cpu;
        if(CParameter.cpu <= 0){
            cpu = Runtime.getRuntime().availableProcessors();
        }else{
            cpu = CParameter.cpu;
        }
        ExecutorService fixedThreadPool = Executors.newFixedThreadPool(cpu);

        /**
         * key: spectrum title, value: an ArrayList of Peptide objects
         */
        ConcurrentHashMap<String,RefPepMatchResult> ref_res = new ConcurrentHashMap<>(spectra2score.size(),0.75F, cpu);

        DBSearchWorker.peptideIndexMap = peptideIndexMap;
        Cloger.getInstance().logger.info("Searching reference database for spectra:"+SpectraInput.spectraMap.size());
        for(String title : SpectraInput.spectraMap.keySet()){
            ref_res.put(title,new RefPepMatchResult());
            ref_res.get(title).score = spectra2score.get(title);
            ref_res.get(title).spectrum_title = title;
            fixedThreadPool.execute(new DBSearchWorker(ref_res.get(title)));
        }

        fixedThreadPool.shutdown();

        try {
            fixedThreadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        // output data
        String ref_detail_file = CParameter.outdir + "/detail.txt";

        BufferedWriter dWriter = new BufferedWriter(new FileWriter(ref_detail_file));
        dWriter.write("spectrum_title\tpeptide\tmodification\texp_mass\tpep_mass\ttol_ppm\ttol_da\tisotope_error\tscore\n");

        for(String spectrum_title:ref_res.keySet()){
            if(ref_res.get(spectrum_title).out_data.size()>=1) {
                try {
                    dWriter.write(StringUtils.join(ref_res.get(spectrum_title).out_data,"\n")+"\n");
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
        }
        dWriter.close();

        if(CParameter.generate_psm_annotation_data) {

            String psm_anno_file = CParameter.outdir + "/psm_annotation.txt";
            BufferedWriter annoWriter = new BufferedWriter(new FileWriter(psm_anno_file));
            annoWriter.write(PeakAnnotation.getPeakAnnotationHeader()+"\n");
            for(String spectrum_title:ref_res.keySet()){
                if(ref_res.get(spectrum_title).psm_annotation_data.size()>=1) {
                    annoWriter.write(StringUtils.join(ref_res.get(spectrum_title).psm_annotation_data,"\n")+"\n");
                }
            }
            annoWriter.close();
        }

        for (PeptideInput peptideInput : peptideInputs) {
            for (JPeptide jPeptide : peptideInput.getPtmIsoforms()) {
                for (int i = 0; i < jPeptide.spectraIndexs.size(); i++) {
                    String spectrum_title = jPeptide.spectraIndexs.get(i);
                    ScoreResult scoreResult = jPeptide.scores.get(i);
                    scoreResult.n_matched_ref_peptides = ref_res.get(spectrum_title).n_matched_peptides;
                    scoreResult.n_matched_ref_peptides_with_better_score = ref_res.get(spectrum_title).n_matched_peptides_with_better_score;
                    if(ref_res.get(spectrum_title).n_matched_peptides_with_better_score>=1){
                        // failed reference peptide filtering
                        jPeptide.valid.add(i,false);
                    }
                }
            }
        }

        long bTime = System.currentTimeMillis();
        Cloger.getInstance().logger.info("Time used for searching reference peptides:" + (bTime - startTime) / 1000 + " s.");

    }


    /**
     * Don't use this function. This is not thread safe. It needs to be improved and fixed.
     */
    public void db_score(String psm_file){

        if(this.ms2peptide.size()>=1) {
            long startTime = System.currentTimeMillis();
            int cpu = 1;
            if (CParameter.cpu <= 0) {
                cpu = Runtime.getRuntime().availableProcessors();
            } else {
                cpu = CParameter.cpu;
            }
            ExecutorService fixedThreadPool = Executors.newFixedThreadPool(cpu);

            /**
             * key: spectrum title, value: an ArrayList of Peptide objects
             */
            // ConcurrentHashMap<String,ArrayList<Peptide>> ms2peptide = new ConcurrentHashMap<>(SpectraInput.spectraMap.size(),0.75F, cpu);
            Long i = 0L;
            for (String title : this.ms2peptide.keySet()) {
                i = i + this.ms2peptide.get(title).size();
                fixedThreadPool.execute(new SinglePSMScoringWorker(this.ms2peptide.get(title), SpectraInput.spectraMap.get(title)));
            }

            fixedThreadPool.shutdown();

            try {
                fixedThreadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            long bTime = System.currentTimeMillis();
            Cloger.getInstance().logger.info("Time used for scoring of reference peptide spectrum matching:" + i + " PSMs, " + (bTime - startTime) / 1000 + " s.");
        }
    }


    public void readDB_v0() throws IOException, InterruptedException {
        File dbFile = new File(this.db);
        FASTAFileReader reader = new FASTAFileReaderImpl(dbFile);
        FASTAElementIterator it = reader.getIterator();

        // variable modifications
        ModificationFactory ptmFactory = ModificationFactory.getInstance();
        /**
         ArrayList<Modification> varPTMs = new ArrayList<>();
         for(String ptmName: varModifications){
         varPTMs.add(ptmFactory.getModification(ptmName));
         }**/

        // digest protein
        //DigestionPreferences digestionPreferences = DigestionPreferences.getDefaultPreferences();
        Enzyme enzyme = getEnzymeByIndex(CParameter.enzyme);
        DigestionParameters digestionPreferences = getDigestionPreferences(enzyme.getName(),CParameter.maxMissedCleavages);

        ArrayList<String> fixM = new ArrayList<>();
        IteratorFactory iteratorModifications = new IteratorFactory(fixM);

        Cloger.getInstance().logger.info("Enzyme: "+enzyme.getName()+", maxMissedCleavages: "+CParameter.maxMissedCleavages);

        // prepare spectrum index for fast peptide query
        HashMap<Integer,ArrayList<String>> spectraIndex = new HashMap<>();
        // step: 0.1
        for(String title : SpectraInput.spectraMap.keySet()){
            Spectrum mSnSpectrum = SpectraInput.spectraMap.get(title);
            double mass = mSnSpectrum.getPrecursor().getMass(mSnSpectrum.getPrecursor().possibleCharges[0]);
            int va = (int) Math.round(10.0*mass);
            if(spectraIndex.containsKey(va)){
                spectraIndex.get(va).add(mSnSpectrum.getSpectrumTitle());
            }else{
                ArrayList<String> spectra_titles = new ArrayList<>();
                spectra_titles.add(mSnSpectrum.getSpectrumTitle());
                spectraIndex.put(va,spectra_titles);
            }
        }

        // if a peptide is searched before, then it will be omitted.
        HashSet<String> searchedPeptides = new HashSet<>();

        int num = 0;
        ExtendedPeptide peptideWithPosition;

        Spectrum spectrum;

        double peptideMass;
        int va;
        double mass;
        double del;

        while (it.hasNext()) {
            FASTAElement el = it.next();
            el.setLineLength(1);
            String headLine[] = el.getHeader().split("\\s+");
            String proID = headLine[0];
            num++;
            // System.out.println(proID);
            String proSeq = el.getSequence().toUpperCase();

            proSeq = proSeq.replaceAll("\\*","");

            //fixedThreadPool.execute(new DBdigestProteinWorker(spectraIndex, proID, proSeq, digestionPreferences, ms2peptide, digestedPeptides));
            //SequenceIterator sequenceIterator = iteratorModifications.getSequenceIterator(proSeq, digestionPreferences, 600.0, 5000.0);
            if(proSeq.length() < CParameter.minPeptideLength){
                continue;
            }
            SequenceIterator sequenceIterator = iteratorModifications.getSequenceIterator(proSeq, digestionPreferences, CParameter.minPeptideMass, CParameter.maxPeptideMass);



            while ((peptideWithPosition = sequenceIterator.getNextPeptide()) != null) {

                Peptide peptide = peptideWithPosition.peptide;


                if (peptide.getSequence().length() < CParameter.minPeptideLength || peptide.getSequence().length() > CParameter.maxPeptideLength) {
                    continue;
                }

                if(searchedPeptides.contains(peptide.getSequence())){
                    continue;
                }else{
                    searchedPeptides.add(peptide.getSequence());
                }

                // Consider modification
                ArrayList<Peptide> peptideIsoforms = PeptideInput.calcPeptideIsoforms(peptide.getSequence());
                //peptideIsoforms.add(peptide);
                for(Peptide p: peptideIsoforms){
                    JPeptide jPeptide = new JPeptide(p);
                    peptideMass = jPeptide.getMass();
                    va = (int) Math.round(10.0*peptideMass);
                    // +/- 0.1 Da
                    for(int i=(va-1);i<=(va+1);i++){
                        if(spectraIndex.containsKey(i)){
                            for(String title : spectraIndex.get(i)){
                                spectrum = SpectraInput.spectraMap.get(title);
                                mass = spectrum.getPrecursor().getMass(spectrum.getPrecursor().possibleCharges[0]);
                                del = Math.abs(peptideMass - mass);
                                if (CParameter.tolu.equalsIgnoreCase("ppm")) {
                                    del = (1.0e6) * 1.0 * del / peptideMass;
                                }
                                if (del <= CParameter.tol) {
                                    //System.out.println("spectra match:"+title+"\t"+mass+"\t"+peptideMass+"\t"+p.getSequence()+"\t"+ ModificationDB.getModificationString(p));
                                    if(ms2peptide.containsKey(title)){
                                        ms2peptide.get(title).add(new JPeptide(p));
                                    }else{
                                        ArrayList<JPeptide> peps = new ArrayList<>();
                                        peps.add(new JPeptide(p));
                                        ms2peptide.put(title,peps);
                                    }
                                }

                            }
                        }
                    }
                }

            }


            if( (num%1000)==0){
                //Cloger.getInstance().logger.info("Finished proteins:"+num);

            }


        }
        reader.close();
        spectraIndex = null;
        Cloger.getInstance().logger.info("Protein sequences:"+num);

    }


    public void readTag() throws IOException {

        /**
        // variable modifications
        ModificationFactory ptmFactory = ModificationFactory.getInstance();
        ArrayList<Modification> varPTMs = new ArrayList<>();
        for(String ptmName: varModifications){
            varPTMs.add(ptmFactory.getModification(ptmName));
        }
         **/


        // Get candidate peptides from tag file
        HashSet<String> spectraIDs = new HashSet<>(SpectraInput.spectraMap.keySet());

        BufferedReader tReader = null;
        try {
            tReader = new BufferedReader(new FileReader(new File(tagFile)));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        String line;
        int n_tags = 0;

        HashMap<String, ArrayList<String>> spectrum2peptides = new HashMap<>();

        while((line = tReader.readLine())!=null){
            String d[] = line.trim().split("\t");

            if(spectraIDs.contains(d[0])){

                String tags[] = d[1].split(";");
                ArrayList<String> tagList = new ArrayList<>();
                for(String tag : tags){
                    tagList.add(tag);
                }
                spectrum2peptides.put(d[0],tagList);
                n_tags = n_tags + tagList.size();


            }

        }

        tReader.close();


        Spectrum spectrum;

        double peptideMass;
        // int va;
        double mass;
        double del;

        for(String title: SpectraInput.spectraMap.keySet()){

            if(!spectrum2peptides.containsKey(title)){
                continue;
            }
            ArrayList<String> tags = spectrum2peptides.get(title);
            Cloger.getInstance().logger.info("tag"+"\t"+title+"\t"+ StringUtils.join(tags,";"));
            for(String tag:tags){
                /**
                // build d peptide object
                Peptide peptide = new Peptide(tag);
                PeptideInput peptideInput = new PeptideInput(peptide);
                peptideInput.addFixedModification(CParameter.fixMods);
                peptideMass = peptideInput.jPeptide.getMass();
                // va = (int) Math.round(10.0*peptideMass);

                spectrum = SpectraInput.spectraMap.get(title);
                mass = spectrum.getPrecursor().getMass(spectrum.getPrecursor().getPossibleCharges().get(0));
                del = Math.abs(peptideMass - mass);
                if (CParameter.tolu) {
                    del = (1.0e6) * 1.0 * del / peptideMass;
                }
                System.out.println(mass+"\t"+peptideMass+"\t"+del);
                if (del <= CParameter.tol) {
                    System.out.println("Yes:"+tag);
                    if(ms2peptide.containsKey(title)){
                        ms2peptide.get(title).add(peptide);
                    }else{
                        ArrayList<Peptide> peps = new ArrayList<>();
                        peps.add(peptide);
                        ms2peptide.put(title,peps);
                    }
                }**/

                // Consider modification
                //ArrayList<Peptide> peptideIsoforms = PeptideInput.calcPeptideIsoforms(peptide,varPTMs,CParameter.maxVarMods);
                ArrayList<Peptide> peptideIsoforms = PeptideInput.calcPeptideIsoforms(tag);
                //peptideIsoforms.add(peptide);
                for(Peptide p: peptideIsoforms){
                    peptideMass = JPeptide.getMass(p);
                    spectrum = SpectraInput.spectraMap.get(title);
                    mass = spectrum.getPrecursor().getMass(spectrum.getPrecursor().possibleCharges[0]);
                    del = Math.abs(peptideMass - mass);
                    if (CParameter.tolu.equalsIgnoreCase("ppm")) {
                        del = (1.0e6) * 1.0 * del / peptideMass;
                    }
                    if (del <= CParameter.tol) {
                        if(ms2peptide.containsKey(title)){
                            ms2peptide.get(title).add(new JPeptide(p));
                        }else{
                            ArrayList<JPeptide> peps = new ArrayList<>();
                            peps.add(new JPeptide(p));
                            ms2peptide.put(title,peps);
                        }
                    }
                }


            }
        }

    }


    /**
     * Calculate the left and right values for specified tol.
     * @param msMass Mass of MS/MS spectrum
     * @param tol tol
     * @param isPpm true if the unit of tol is ppm
     * @return
     */
    public double[] getRangeOfMass(double msMass, double tol, boolean isPpm){
        double [] massRange = new double [2];
        if(isPpm){
            massRange[0] = 1.0*msMass/(1+1.0*tol/(1.0e6));
            massRange[1] = 1.0*msMass/(1-1.0*tol/(1.0e6));

        }else{
            massRange[0] = msMass-tol;
            massRange[1] = msMass+tol;
        }
        return massRange;
    }


    public ArrayList<Peptide> getPeptideFromSQL(ResultSet rs) throws SQLException, IOException, ClassNotFoundException {
        ArrayList<Peptide> peptides = new ArrayList<>();
        // peptide sequence + modification
        HashSet<String> pepmod = new HashSet<>();
        String pmod;
        while(rs.next()){
            pmod = rs.getString("pep_mod") + ":" + rs.getString("modification");
            if(pepmod.contains(pmod)){
                continue;
            }else{
                pepmod.add(pmod);
            }
            byte buf[] = rs.getBytes("peptide");
            ObjectInputStream objectInputStream = new ObjectInputStream(new ByteArrayInputStream(buf));
            Peptide p = (Peptide) objectInputStream.readObject();
            peptides.add(p);

        }

        return peptides;

    }


    public void readDBSQL() throws IOException, SQLException, ClassNotFoundException {

        // read digested peptides from SQL database
        Connection connection = DriverManager.getConnection("jdbc:sqlite:" + this.db);
        PreparedStatement pstmt = connection.prepareStatement("select sequence,peptide,pep_mod,modification from prodb where mass >= ? and mass <= ?");
        String psql;
        int npep = 0;
        Spectrum spectrum;
        for(String title : SpectraInput.spectraMap.keySet()){
            spectrum = SpectraInput.spectraMap.get(title);
            double mass = spectrum.getPrecursor().getMass(spectrum.getPrecursor().possibleCharges[0]);
            double massRange[] = getRangeOfMass(mass,CParameter.tol,CParameter.tolu.equalsIgnoreCase("ppm"));
            pstmt.setDouble(1,massRange[0]);
            pstmt.setDouble(2,massRange[1]);
            ResultSet rs    = pstmt.executeQuery();

            ArrayList<Peptide> rsPeptide = getPeptideFromSQL(rs);
            npep = npep + rsPeptide.size();
            /**
            for(Peptide pp:rsPeptide) {
                System.out.println("spectra match:" + title + "\t" + mass + "\t"+JPeptide.getMass(pp)+"\t" + pp.getSequence() + "\t" + ModificationDB.getModificationString(pp));
            }**/

            ArrayList<JPeptide> jp = new ArrayList<>();
            for(Peptide p : rsPeptide){
                jp.add(new JPeptide(p));
            }
            ms2peptide.put(spectrum.getSpectrumTitle(),jp);
        }
        pstmt.close();
        connection.close();
        Cloger.getInstance().logger.info("Peptides from reference database:"+npep);



    }


    public Comparator<Spectrum> comparator = (s1, s2) -> {
        if (s2.getPrecursor().getMass(s2.getPrecursor().possibleCharges[0]) >
                s1.getPrecursor().getMass(s1.getPrecursor().possibleCharges[0])) {
            return -1;
        } else if (s2.getPrecursor().getMass(s2.getPrecursor().possibleCharges[0]) ==
                s1.getPrecursor().getMass(s1.getPrecursor().possibleCharges[0])) {
            return 0;
        } else {
            return 1;
        }

    };





}
