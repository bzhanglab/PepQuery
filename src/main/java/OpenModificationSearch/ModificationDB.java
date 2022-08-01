package main.java.OpenModificationSearch;

import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.experiment.identification.protein_sequences.SingleProteinSequenceProvider;
import com.compomics.util.experiment.identification.utils.ModificationUtils;
import com.compomics.util.experiment.io.biology.protein.SequenceProvider;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import com.compomics.util.experiment.biology.modifications.*;
import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import com.compomics.util.pride.CvTerm;
import main.java.pg.*;
import main.java.PSMMatch.JPeptide;
import main.java.PSMMatch.ScoreResult;
import main.java.util.Cloger;
import org.apache.commons.lang3.StringUtils;
import org.paukov.combinatorics3.Generator;
import uk.ac.ebi.pride.utilities.pridemod.io.unimod.model.Specificity;
import uk.ac.ebi.pride.utilities.pridemod.io.unimod.model.Unimod;
import uk.ac.ebi.pride.utilities.pridemod.io.unimod.model.UnimodModification;
import uk.ac.ebi.pride.utilities.pridemod.io.unimod.xml.UnimodReader;
import javax.xml.bind.JAXBException;
import java.io.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;



/**
 * 1. Build peptide database including peptide sequences and fixed modifications
 * 2. Build index for modifications
 *
 * Query method:
 * Define a windows (250 Da), get all the peptides with mass error less than the defined window compared with query
 * spectrum. Only consider the peptides with mass error less than pre-defined mass error (for example 10 ppm). For this
 * search, currently consider one variable modification for each PSM.
 */
public class ModificationDB {

    private static ModificationDB instance = null;
    private ModificationFactory ptmFactory = null;

    public static boolean save_mod2file = false;

    private ModificationDB(){
        ptmFactory = ModificationFactory.getInstance();
        try {
            importPTMsFromUnimod();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static ModificationDB getInstance() {
        if (instance == null) {
            instance = new ModificationDB();
        }
        return instance;
    }

    /**
     * Peptide index for spectrum matching
     * The peptides are already considered possible fixed and variable modifications
     */
    private HashMap<Long, ArrayList<JPeptide>> peptideIndexMap;

    public void add_peptide_index(HashMap<Long, ArrayList<JPeptide>> peptide_index){
        this.peptideIndexMap = peptide_index;
    }
    public HashMap<Long, ArrayList<JPeptide>> get_peptide_index(){
        return this.peptideIndexMap;
    }

    public void clear_peptide_index(){
        if(this.get_peptide_index()!=null) {
            this.peptideIndexMap.clear();
        }
    }

    private boolean debug = false;

    /**
     * Integer value of a mass
     */
    private HashMap<Integer, Integer> mass2index = new HashMap<>();
    private HashMap<Integer, Integer> mass2indexUp = new HashMap<>();



    // private static SequenceMatchingPreferences sequenceMatchingPreferences = new SequenceMatchingPreferences();
    // private static SequenceMatchingPreferences ptmSequenceMatchingPreferences = new SequenceMatchingPreferences();


    private HashMap<String, String> titleTFullName = new HashMap<>();

    /**
     * Store PTMs that will be considered in peptide matching
     */
    private  ArrayList<Modification> ptms = new ArrayList<>();
    private  boolean isPTMsort = false;


    /**
     * Index for modification searching
     * key=Integer value of a mass, value=all modifications to this key
     * Only the modifications from UniMod database.
     */
    private HashMap<Long,ArrayList<Modification>> ptmIndex = new HashMap<>();


    private int offset = 200;

    public static double leftMass = -250.0;
    public static double rightMass = 250.0;

    //
    public static String out_dir = "./";


    public static void main(String[] args) {

    }

    /**
     * Get Modifications whose masses are matched with delta mass (mass).
     * @param ptm_mass delta mass
     * @param precursorMass Peptide mass not MS/MS spectrum precursor mass. This is the peptide mass before occurring the variable modification.
     * @return An array of Modification objects.
     */
    private ArrayList<Modification> getCandidatePTMs(double ptm_mass, double precursorMass){

        ArrayList<Modification> result = new ArrayList<>();

        long ptmmass = Math.round(ptm_mass);
        if(ptmIndex.containsKey(ptmmass) || ptmIndex.containsKey(ptmmass-1) || ptmIndex.containsKey(ptmmass+1)) {
            // left: 1Da, right: 1Da
            for(long i=ptmmass-1;i<=ptmmass+1;i++){
                if(ptmIndex.containsKey(i)){
                    for(Modification ptm:ptmIndex.get(i)){
                        double mdel = ptm_mass - ptm.getMass();
                        if (CParameter.tolu.equalsIgnoreCase("ppm")) {
                            // precursorMass+ptm.getMass() = final peptide mass
                            // please note this is not theoretical peptide mass
                            // so the real tol could exceed CParameter.tol even for now
                            // (Math.abs(mdel) <= CParameter.tol) is true
                            // so after using getCandidatePTMs, we need to
                            // check this condition again.
                            mdel = mdel / (precursorMass+ptm.getMass()) * 1000000.0;
                        }

                        if (Math.abs(mdel) <= CParameter.tol) {
                            result.add(ptm);
                        }
                    }
                }
            }
        }
        return result;

    }

    /**
     * Get all Peptides matched to the spectrum. Find if the peptide could be matched to the target spectrum with any
     * modifications from UniMod.
     * @param spectrum A Spectrum object
     * @param peptide A Peptide object. The peptide is already considered fixed modification if applicable.
     * @return An array of Peptide objects
     */
    public ArrayList<Peptide> getCandidatePTMs(Spectrum spectrum, Peptide peptide){

        ArrayList<Peptide> modPeptides = new ArrayList<>();


        HashSet<Integer> fixedModSites = new HashSet<>();
        ModificationMatch modificationMatchs[] = peptide.getVariableModifications();
        for(ModificationMatch match: modificationMatchs){
            fixedModSites.add(match.getSite());
            // debug
            //System.out.println("PTM site:"+match.getSite()+" -> "+match.getModification());
        }

        double precursorMass = spectrum.getPrecursor().getMass(spectrum.getPrecursor().possibleCharges[0]);
        // isotope_masses store mono mass for precursor
        ArrayList<Double> isotope_masses = CParameter.get_precursor_mass_with_isotope_error(precursorMass);

        double pep_mass = JPeptide.getMass(peptide);
        for(double pre_mass:isotope_masses) {
            double delmass = pre_mass - pep_mass;
            //ArrayList<Modification> rptms = getCandidatePTMs(delmass, pre_mass);
            ArrayList<Modification> rptms = getCandidatePTMs(delmass, pep_mass);
            Iterator<Modification> iterator = rptms.iterator();
            while (iterator.hasNext()) {
                // for each PTM
                Modification ptm = iterator.next();

                if (ptm.getModificationType() == ModificationType.modn_protein || ptm.getModificationType() == ModificationType.modnaa_protein ||
                        ptm.getModificationType() == ModificationType.modc_protein || ptm.getModificationType() == ModificationType.modcaa_protein) {

                    iterator.remove();
                } else {

                    // possibleSites = peptide.getPotentialModificationSites(ptm, sequenceMatchingPreferences, ptmSequenceMatchingPreferences);
                    // ArrayList<Integer> possibleSites = peptide.getPotentialModificationSitesNoCombination(ptm, peptide.getSequence(),0);
                    int pSites[] = ModificationUtils.getPossibleModificationSites(peptide, ptm, sequenceProvider, modificationsSequenceMatchingParameters);

                    if (pSites.length >= 1) {

                        ArrayList<Integer> possibleSites = new ArrayList<>();
                        for (int i : pSites) {
                            possibleSites.add(i);
                        }
                        possibleSites.removeAll(fixedModSites);

                        if (possibleSites != null && possibleSites.size() >= 1) {


                            HashMap<Integer, Modification> position2ptm = new HashMap<>();

                            for (Integer i : possibleSites) {
                                position2ptm.put(i, ptm);
                            }

                            if (position2ptm.size() >= 1) {
                                // >=1 possible modification site
                                Set<Integer> ptmPosition = position2ptm.keySet();
                                // the max number of modifications to consider
                                Iterator<List<Integer>> iter = Generator.combination(ptmPosition).simple(1).iterator();
                                while (iter.hasNext()) {
                                    List<Integer> iCombination = iter.next();
                                    // Peptide modPeptide = new Peptide(peptide.getSequence(), peptide.getModificationMatches());
                                    // getVariableModifications: here it is fixed modification
                                    Peptide modPeptide = new Peptide(peptide.getSequence(), peptide.getVariableModifications());
                                    for (Integer pos : iCombination) {
                                        // pos : 1-based position
                                        modPeptide.addVariableModification(new ModificationMatch(position2ptm.get(pos).getName(), pos));
                                    }
                                    modPeptides.add(modPeptide);

                                }
                            }
                        }
                    }
                }


            }
        }
        return(modPeptides);

    }


    /**
     *
     * @param peptide
     * @param ptm
     * @return
     */
    public ArrayList<Peptide> getCandidatePTMs(Peptide peptide, Modification ptm){

        ArrayList<Peptide> modPeptides = new ArrayList<>();


        HashSet<Integer> fixedModSites = new HashSet<>();
        ModificationMatch modificationMatchs[] = peptide.getVariableModifications();
        for(ModificationMatch match: modificationMatchs){
            fixedModSites.add(match.getSite());
            // debug
            //System.out.println("PTM site:"+match.getSite()+" -> "+match.getModification());
        }

        if (ptm.getModificationType() == ModificationType.modn_protein || ptm.getModificationType() == ModificationType.modnaa_protein ||
                ptm.getModificationType() == ModificationType.modc_protein || ptm.getModificationType() == ModificationType.modcaa_protein) {

            return modPeptides;
        } else {

            // possibleSites = peptide.getPotentialModificationSites(ptm, sequenceMatchingPreferences, ptmSequenceMatchingPreferences);
            // ArrayList<Integer> possibleSites = peptide.getPotentialModificationSitesNoCombination(ptm, peptide.getSequence(),0);
            int pSites[] = ModificationUtils.getPossibleModificationSites(peptide, ptm, sequenceProvider, modificationsSequenceMatchingParameters);

            if (pSites.length >= 1) {

                ArrayList<Integer> possibleSites = new ArrayList<>();
                for (int i : pSites) {
                    possibleSites.add(i);
                }
                possibleSites.removeAll(fixedModSites);

                if (possibleSites != null && possibleSites.size() >= 1) {


                    HashMap<Integer, Modification> position2ptm = new HashMap<>();

                    for (Integer i : possibleSites) {
                        position2ptm.put(i, ptm);
                    }

                    if (position2ptm.size() >= 1) {
                        // >=1 possible modification site
                        Set<Integer> ptmPosition = position2ptm.keySet();
                        // the max number of modifications to consider
                        Iterator<List<Integer>> iter = Generator.combination(ptmPosition).simple(1).iterator();
                        while (iter.hasNext()) {
                            List<Integer> iCombination = iter.next();
                            // Peptide modPeptide = new Peptide(peptide.getSequence(), peptide.getModificationMatches());
                            // getVariableModifications: here it is fixed modification
                            Peptide modPeptide = new Peptide(peptide.getSequence(), peptide.getVariableModifications());
                            for (Integer pos : iCombination) {
                                // pos : 1-based position
                                modPeptide.addVariableModification(new ModificationMatch(position2ptm.get(pos).getName(), pos));
                            }
                            modPeptides.add(modPeptide);

                        }
                    }
                }
            }
        }
        return(modPeptides);

    }




    /**
     * Import modification information from Unimod.
     */
    private void importPTMsFromUnimod() throws IOException {
        //ptmFactory.clearFactory();
        ArrayList<String > residues;

        ArrayList<String> testSite = new ArrayList<>();
        ArrayList<String> testPosition = new ArrayList<>();

        InputStream inputStream = ModificationDB.class.getResourceAsStream("/main/resources/unimod.xml");
        // main/java/OpenModificationSearch/unimod.xml

        UnimodReader unimodreader = null;
        try {
            unimodreader = new UnimodReader(inputStream);
        } catch (JAXBException e) {
            e.printStackTrace();
        }

        Unimod unimod = unimodreader.getUnimodObject();

        List<UnimodModification> unimodModificationList = unimod.getModifications().getMod();

        Cloger.getInstance().logger.info("All modifications in unimod:"+unimodModificationList.size());

        HashMap<String,Integer> modificationClassMap = new HashMap<>();


        ArrayList<String> mod2file_data = new ArrayList<>();

        //Parsing and get each modification;
        int i=0;
        for (UnimodModification unimodModification: unimodModificationList){

            Double monoMass = unimodModification.getDelta().getMonoMass().doubleValue();//Get modification Mass

            if(monoMass < leftMass || monoMass > rightMass){
                continue;
            }

            String modificationTitle = unimodModification.getTitle();
            String modificationName = unimodModification.getFullName();
            titleTFullName.put(modificationTitle, modificationName);

            String unimod_accession = "UNIMOD:"+String.valueOf(unimodModification.getRecordId());

            //Get amino acid modification detailed modification;
            List<Specificity> specificityList = unimodModification.getSpecificity();
            //String fullName = unimodModification.getFullName();

            for(Specificity specificity: specificityList){
                String site = specificity.getSite();
                String position = specificity.getPosition();
                String classification = specificity.getClassification();

                // Only importing these modifications when performing target protein (known) identification. But it's optional.
                // There is a global parameter in CParameter "addAAsubstitutionMods" to control this.
                // Only when CParameter.addAAsubstitutionMods is false, then it will be ignored.
                if(classification.equalsIgnoreCase("AA substitution") && !CParameter.addAAsubstitutionMods){
                    continue;
                }
                //System.out.println(classification);

                if(modificationClassMap.containsKey(classification)){
                    modificationClassMap.put(classification,modificationClassMap.get(classification)+1);
                }else{
                    modificationClassMap.put(classification,1);
                }

                // output
                if(save_mod2file){
                    mod2file_data.add(modificationTitle + " of " + site+"\t"+modificationName+"\t"+monoMass+"\t"+site+"\t"+position+"\t"+classification);
                }


                Modification ptm = null;
                String ptmName;
                i++;

                if(site.equals("N-term")){

                    if(position.equals("Any N-term")){
                        ptmName = modificationTitle + " of " + site;
                        //ptm = new Modification(ModificationType.modn_peptide, ptmName, monoMass, null);
                        ptm = new Modification(ModificationType.modn_peptide,ptmName,monoMass,null,ModificationCategory.Other);
                        ptm.setShortName(modificationTitle);
                    }else if(position.equals("Protein N-term")){
                        ptmName = modificationTitle + " of protein " + site;
                        ptm = new Modification(ModificationType.modn_protein,ptmName,monoMass,null,ModificationCategory.Other);
                        ptm.setShortName(modificationTitle);
                    }else{

                        if(debug) {

                            Cloger.getInstance().logger.info("position:" + position);
                            Cloger.getInstance().logger.info("mass:" + monoMass);
                            Cloger.getInstance().logger.info("site:" + site);
                        }
                    }
                }else if(site.equals("C-term")){
                    if(position.equals("Any C-term")){
                        ptmName = modificationTitle + " of " + site;
                        ptm = new Modification(ModificationType.modc_peptide, ptmName,monoMass,null,ModificationCategory.Other);
                        ptm.setShortName(modificationTitle);
                    }else if(position.equals("Protein C-term")){
                        ptmName = modificationTitle + " of protein " + site;
                        ptm = new Modification(ModificationType.modc_protein, ptmName,monoMass,null,ModificationCategory.Other);
                        ptm.setShortName(modificationTitle);
                    }
                }else {
                    residues = new ArrayList<>();
                    residues.add(site);
                    if(position.equals("Any N-term")){
                        ptmName = modificationTitle + " of " + site;
                        // Modification at the N terminus of a peptide at particular amino acids.
                        ptm = new Modification(ModificationType.modnaa_peptide, ptmName, monoMass, residues,ModificationCategory.Other);
                        ptm.setShortName(modificationTitle);
                    }else if(position.equals("Protein N-term")){
                        ptmName = modificationTitle + " of protein " + site;
                        // Modification at the N terminus of a protein at particular amino acids.
                        ptm = new Modification(ModificationType.modnaa_protein, ptmName, monoMass, residues,ModificationCategory.Other);
                        ptm.setShortName(modificationTitle);
                    }else if (position.equals("Any C-term")){
                        ptmName = modificationTitle + " of " + site;
                        ptm = new Modification(ModificationType.modcaa_peptide, ptmName, monoMass, residues,ModificationCategory.Other);
                        ptm.setShortName(modificationTitle);
                    }else if(position.equals("Protein C-term")){
                        ptmName = modificationTitle + " of protein " + site;
                        ptm = new Modification(ModificationType.modcaa_protein, ptmName, monoMass, residues,ModificationCategory.Other);
                        ptm.setShortName(modificationTitle);
                    }else {
                        ptmName = modificationTitle + " of " + site;
                        ptm = new Modification(ModificationType.modaa, ptmName, monoMass, residues,ModificationCategory.Other);
                        ptm.setShortName(modificationTitle);
                    }
                }

                //System.out.println("The ptm name is "+ptmName+" and mass is "+ptm.getMass());

                if(ptm!=null) {

                    // This is important. Indicating this modification is a amino acid substitution.
                    if(classification.equalsIgnoreCase("AA substitution") && CParameter.addAAsubstitutionMods){
                        ptm.setCategory(ModificationCategory.Nucleotide_Substitution_One);
                    }
                    CvTerm cvTerm = new CvTerm();
                    cvTerm.setAccession(unimod_accession);
                    ptm.setUnimodCvTerm(cvTerm);
                    if(ptmFactory.containsModification(ptm.getName())){
                        // if the ptm exists in default PTM list, don't use the ptm information from UniMod
                        ptms.add(ptmFactory.getModification(ptm.getName()));
                    }else{
                        ptms.add(ptm);
                        ptmFactory.addUserModification(ptm);
                    }


                    // Build index for modification search
                    long ptmmass = Math.round(monoMass);
                    if(ptmIndex.containsKey(ptmmass)){
                        ptmIndex.get(ptmmass).add(ptm);
                    }else{
                        ArrayList<Modification> pList = new ArrayList<>();
                        pList.add(ptm);
                        ptmIndex.put(ptmmass,pList);
                    }
                }
                //System.out.println(ptm.getMass());

                if(!testSite.contains(site)){
                    testSite.add(site);
                }
                if(!testPosition.contains(position)){
                    testPosition.add(position);
                }

            }
        }

        if(save_mod2file){
            // Export detailed modification information
            String out_mod_file = out_dir + "/Unimod.tsv";
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(out_mod_file)));
            bw.write("mod_title\tmod_name\tmono_mass\tsite\tposition\tclassification\n");
            for(String line: mod2file_data) {
                bw.write(line+"\n");
            }
            bw.close();
        }

        if(debug) {

            int total_mod = 0;
            for (String classification : modificationClassMap.keySet()) {
                Cloger.getInstance().logger.info(classification + "\t" + modificationClassMap.get(classification));
                total_mod = total_mod + modificationClassMap.get(classification);
            }

            Cloger.getInstance().logger.info("Total modifications imported:" + this.ptms.size());
        }

    }


    /**
     *
     */
    private Comparator<Modification> comparator = new Comparator<Modification>() {

        public int compare(Modification s1, Modification s2) {
            if (s1.getMass() < s2.getMass()) {
                return -1;
            } else if (s1.getMass() > s2.getMass()) {
                return 1;
            } else {
                return 0;
            }
        }

    };


    private synchronized double getMinModificationMass(){
        double mass = 0.0;
        if(this.ptms.size()>=1){
            if(this.isPTMsort){
                mass = this.ptms.get(0).getMass();
            }else{
                Collections.sort(this.ptms, comparator);
                this.isPTMsort = true;
                mass = this.ptms.get(0).getMass();
            }
        }else{
            System.err.println("There is not modification!");
            System.exit(0);
        }
        return(mass);
    }




    /**
     * Digest proteins with considering fixed modifications and then build peptide index for open search
     * @param db database file in FASTA format
     * @param fixedModifications Fixed modification
     * @return HashMap<Long,PeptideIndex>

    private HashMap<Long,PeptideIndex> buildPeptideDB(String db, ArrayList<String> fixedModifications, HashSet<String> target_peptides){

        long startTime=System.currentTimeMillis();

        HashMap<Long,PeptideIndex> peptideIndexMap = new HashMap<>();

        String pepIndex_file = db+ ".pepIndex";
        File PEPInd = new File(pepIndex_file);
        boolean useIndex = false;
        if(PEPInd.isFile() && useIndex){

            Cloger.getInstance().logger.info("Read peptide index file: "+pepIndex_file);

            // read index object
            try {
                FileInputStream pepindex_fi = new FileInputStream(new File(pepIndex_file));
                ObjectInputStream pepindex_oi = new ObjectInputStream(pepindex_fi);
                peptideIndexMap = (HashMap<Long,PeptideIndex>) pepindex_oi.readObject();

            } catch (FileNotFoundException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            } catch (ClassNotFoundException e) {
                e.printStackTrace();
            }
        }else {

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
            DigestionParameters digestionParameters = DatabaseInput.getDigestionPreferences(enzyme.getName(),CParameter.maxMissedCleavages);

            //DigestionPreferences digestionPreferences = DigestionPreferences.getDefaultPreferences();
            // if a peptide is searched before, then it will be omitted.
            HashSet<String> searchedPeptides = new HashSet<>();


            int cpu = Runtime.getRuntime().availableProcessors();
            ExecutorService fixedThreadPool = Executors.newFixedThreadPool(cpu);

            //ConcurrentHashMap<String, HashMap<String, Double>> pro2pep = new ConcurrentHashMap<>();
            ConcurrentHashMap<String, ArrayList<Peptide>> pro2pep = new ConcurrentHashMap<>();

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
                    fixedThreadPool.execute(new DigestProteinWorker(proID, proSeq, fixedModifications, digestionParameters, pro2pep));
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
            Cloger.getInstance().logger.info("Protein sequences:" + num);


            for (String proID : pro2pep.keySet()) {
                for (String pep : pro2pep.get(proID).keySet()) { //at main.java.OpenModificationSearch.ModificationDB.buildPeptideDB(ModificationDB.java:874) Exception in thread "main" java.lang.NullPointerException
                    if (!searchedPeptides.contains(pep)) {
                        if(target_peptides.contains(pep)){
                            continue;
                        }
                        searchedPeptides.add(pep);
                        double mass = pro2pep.get(proID).get(pep);
                        long intValue = Math.round(mass);
                        if (peptideIndexMap.containsKey(intValue)) {
                            peptideIndexMap.get(intValue).peptideSequences.add(pep);
                            peptideIndexMap.get(intValue).peptideMasses.add(mass);
                        } else {
                            PeptideIndex peptideIndex = new PeptideIndex(pep, mass);
                            peptideIndexMap.put(intValue, peptideIndex);

                        }
                    }
                }
            }

            long bTime = System.currentTimeMillis();
            Cloger.getInstance().logger.info("Total unique peptides: "+searchedPeptides.size());
            Cloger.getInstance().logger.info("Time used for building peptide index:" + (bTime - startTime) / 1000 + " s.");


        }

        long endTime = System.currentTimeMillis();

        Cloger.getInstance().logger.info("Time used for peptide indexing:" + (endTime - startTime) / 1000 + " s.");

        return(peptideIndexMap);
    }**/


    /**
    private HashMap<Long,PeptideIndex> buildPeptideDB(String db, HashSet<String> target_peptides){

        long startTime=System.currentTimeMillis();

        HashMap<Long,PeptideIndex> peptideIndexMap = new HashMap<>();

        String pepIndex_file = db+ ".pepIndex";
        File PEPInd = new File(pepIndex_file);
        boolean useIndex = false;
        if(PEPInd.isFile() && useIndex){

            Cloger.getInstance().logger.info("Read peptide index file: "+pepIndex_file);

            // read index object
            try {
                FileInputStream pepindex_fi = new FileInputStream(new File(pepIndex_file));
                ObjectInputStream pepindex_oi = new ObjectInputStream(pepindex_fi);
                peptideIndexMap = (HashMap<Long,PeptideIndex>) pepindex_oi.readObject();

            } catch (FileNotFoundException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            } catch (ClassNotFoundException e) {
                e.printStackTrace();
            }
        }else {

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
            DigestionParameters digestionParameters = DatabaseInput.getDigestionPreferences(enzyme.getName(),CParameter.maxMissedCleavages);

            // DigestionPreferences digestionPreferences = DigestionPreferences.getDefaultPreferences();
            // if a peptide is searched before, then it will be omitted.
            ConcurrentHashMap<String, ArrayList<JPeptide>> searchedPeptides = new ConcurrentHashMap<>();

            // protein digestion
            int cpu = Runtime.getRuntime().availableProcessors();
            ExecutorService fixedThreadPool = Executors.newFixedThreadPool(cpu);

            // ConcurrentHashMap<String, HashMap<String, Double>> pro2pep = new ConcurrentHashMap<>();
            // ConcurrentHashMap<String, ArrayList<Peptide>> pro2pep = new ConcurrentHashMap<>();
            ConcurrentHashMap<String, HashSet<String>> pro2pepSeq = new ConcurrentHashMap<>();

            int num = 0;
            try {
                while (it.hasNext()) {
                    FASTAElement el = it.next();
                    el.setLineLength(1);
                    String headLine[] = el.getHeader().split("\\s+");
                    String proID = headLine[0];
                    num++;
                    String proSeq = el.getSequence();
                    // ArrayList<Peptide> pList = new ArrayList<>(20);
                    // pro2pep.put(proID,pList);
                    pro2pepSeq.put(proID, new HashSet<String>());
                    // digest
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

            for (String proID : pro2pepSeq.keySet()) {
                for (String pep : pro2pepSeq.get(proID)){
                    ArrayList<Peptide> pList = new ArrayList<>();
                    searchedPeptides.put(pep, pList);
                }
            }
            pro2pepSeq.clear();

            Cloger.getInstance().logger.info("Protein sequences:" + num);
            // protein digestion done.

            // generate peptide forms for all unique peptide sequences.
            fixedThreadPool = Executors.newFixedThreadPool(cpu);
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
            // generate peptide forms for all unique peptide sequences finished.

            // generate peptide indexing for fast unrestricted modified peptide matching
            long n_peptide_forms = 0;

            for (String peptide_seq : searchedPeptides.keySet()) {
                for (Peptide pep : searchedPeptides.get(peptide_seq)) { //at main.java.OpenModificationSearch.ModificationDB.buildPeptideDB(ModificationDB.java:874) Exception in thread "main" java.lang.NullPointerException
                    n_peptide_forms = n_peptide_forms + 1;
                    double mass = JPeptide.getMass(pep);
                    long intValue = Math.round(mass);
                    if (peptideIndexMap.containsKey(intValue)) {
                        //peptideIndexMap.get(intValue).peptideSequences.add(pep);
                        peptideIndexMap.get(intValue).peptides.add(pep);
                        // peptideIndexMap.get(intValue).peptideMasses.add(mass);
                    } else {
                        // PeptideIndex peptideIndex = new PeptideIndex(pep, mass);
                        PeptideIndex peptideIndex = new PeptideIndex(pep);
                        peptideIndexMap.put(intValue, peptideIndex);

                    }
                }
            }

            long bTime = System.currentTimeMillis();
            Cloger.getInstance().logger.info("Total unique peptide: "+searchedPeptides.size() + ", total peptide forms:"+n_peptide_forms);
            Cloger.getInstance().logger.info("Time used for building peptide index:" + (bTime - startTime) / 1000 + " s.");

            searchedPeptides.clear();

        }

        long endTime = System.currentTimeMillis();

        Cloger.getInstance().logger.info("Time used for peptide indexing:" + (endTime - startTime) / 1000 + " s.");

        return(peptideIndexMap);
    }**/


    /**
     * Build peptide index for open search. Read candidates from a peptide file.
     * @param tagFile
     * @param spectraIDs
     * @param tagType
     * @return HashMap<String, HashMap<Long, PeptideIndex>>
     */
    private HashMap<String, HashMap<Long, PeptideIndex>> buildPeptideDB(String tagFile, HashSet<String> spectraIDs, int tagType) throws IOException {

        long startTime=System.currentTimeMillis();

        HashMap<String, HashMap<Long, PeptideIndex>> pepIndexMap = new HashMap<>();

        BufferedReader tReader = null;
        try {
            tReader = new BufferedReader(new FileReader(new File(tagFile)));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        String line;
        int n_tags = 0;
        while((line = tReader.readLine())!=null){
            String d[] = line.trim().split("\t");

            if(spectraIDs.contains(d[0])){
                HashMap<Long,PeptideIndex> peptideIndexMap = new HashMap<>();
                String tags[] = d[1].split(";");
                n_tags = n_tags + tags.length;
                for(String tag:tags){
                    // build d peptide object
                    Peptide peptide = new Peptide(tag);
                    PeptideInput.addFixedModification(CParameter.fixMods,peptide);
                    double mass = JPeptide.getMass(peptide);
                    long intValue = Math.round(mass);
                    if (peptideIndexMap.containsKey(intValue)) {
                        peptideIndexMap.get(intValue).peptideSequences.add(tag);
                        peptideIndexMap.get(intValue).peptideMasses.add(mass);
                    } else {
                        PeptideIndex peptideIndex = new PeptideIndex(tag, mass);
                        peptideIndexMap.put(intValue, peptideIndex);

                    }


                }

                pepIndexMap.put(d[0], peptideIndexMap);
            }

        }


        tReader.close();


        long bTime = System.currentTimeMillis();
        double mean_tags = 1.0*n_tags/spectraIDs.size();
        Cloger.getInstance().logger.info("Total unique peptides per spectrum: "+mean_tags);
        Cloger.getInstance().logger.info("Time used for building peptide index:" + (bTime - startTime) / 1000 + " s.");

        return(pepIndexMap);
    }


    /**
     * For each target spectrum, find matched reference peptides with any modifications from UniMod.
     * @param spectrum A Spectrum object
     * @param peptideIndexMap Peptide indexing hash map which contains all reference peptides. The key is peptide mass index (integer value). The value for each key is a PeptideIndex object.
     * @param targetPeptideSeq Target peptide sequence matched to the target spectrum.
     * @return A hash map which contains scoring result.
     */
    private ConcurrentHashMap<Long, ArrayList<ScoreResult>> getCandidateRefPeptides(Spectrum spectrum,
                                                                                    HashMap<Long,ArrayList<JPeptide>> peptideIndexMap,
                                                                                    String targetPeptideSeq) {

        if(leftMass < getMinModificationMass()){
            leftMass = getMinModificationMass();
        }

        if(debug) {
            Cloger.getInstance().logger.info("left mass:" + leftMass);
            Cloger.getInstance().logger.info("left mass:" + rightMass);
        }
        double precursorMass = spectrum.getPrecursor().getMass(spectrum.getPrecursor().possibleCharges[0]);
        if(debug) {
            Cloger.getInstance().logger.info("precursorMass=" + precursorMass + ",leftMass=" + leftMass + ",rightMass=" + rightMass);
        }

        // HyperscoreMatch.generateFactorialValues(60);


        // peptide mass with modifications
        // The max and min delta mass (PTM mass) will be considered.
        long maxIntValue = Math.round(precursorMass-leftMass+2);
        long minIntValue = Math.round(precursorMass-rightMass-2);
        //System.out.println("left mass:" + leftMass);
        //System.out.println("left mass:" + rightMass);
        //System.out.println("precursorMass=" + precursorMass + ",leftMass=" + leftMass + ",rightMass=" + rightMass);

        int cpu = CParameter.cpu;

        if (cpu == 0) {
            cpu = Runtime.getRuntime().availableProcessors();
        }
        //System.out.println("CPU:" + cpu);
        ExecutorService fixedThreadPool = Executors.newFixedThreadPool(cpu);

        ConcurrentHashMap<Long,ArrayList<ScoreResult>> res = new ConcurrentHashMap<>(500,0.75F,cpu);

        //int ind = 1;

        for (long i = minIntValue; i <= maxIntValue; i++) {
            ArrayList<ScoreResult> sList = new ArrayList<>();
            res.put(i,sList);
            //
            fixedThreadPool.execute(new ModPepQueryAndScoringWorker(i,spectrum, peptideIndexMap,res,targetPeptideSeq));
        }

        fixedThreadPool.shutdown();

        try {
            fixedThreadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }


        return(res);

    }

    private ConcurrentHashMap<Long, ArrayList<ScoreResult>> getCandidateRefPeptidesFast(Spectrum spectrum,
                                                                                    HashMap<Long,ArrayList<JPeptide>> peptideIndexMap,
                                                                                    String targetPeptideSeq) {

        if(leftMass < getMinModificationMass()){
            leftMass = getMinModificationMass();
        }

        if(debug) {
            Cloger.getInstance().logger.info("left mass:" + leftMass);
            Cloger.getInstance().logger.info("left mass:" + rightMass);
        }
        double precursorMass = spectrum.getPrecursor().getMass(spectrum.getPrecursor().possibleCharges[0]);
        if(debug) {
            Cloger.getInstance().logger.info("precursorMass=" + precursorMass + ",leftMass=" + leftMass + ",rightMass=" + rightMass);
        }

        // HyperscoreMatch.generateFactorialValues(60);


        // peptide mass with modifications
        // The max and min delta mass (PTM mass) will be considered.
        long maxIntValue = Math.round(precursorMass-leftMass+2);
        long minIntValue = Math.round(precursorMass-rightMass-2);
        //System.out.println("left mass:" + leftMass);
        //System.out.println("left mass:" + rightMass);
        //System.out.println("precursorMass=" + precursorMass + ",leftMass=" + leftMass + ",rightMass=" + rightMass);

        int cpu = CParameter.cpu;

        if (cpu == 0) {
            cpu = Runtime.getRuntime().availableProcessors();
        }
        //System.out.println("CPU:" + cpu);
        ExecutorService fixedThreadPool = Executors.newFixedThreadPool(cpu);

        ConcurrentHashMap<Long,ArrayList<ScoreResult>> res = new ConcurrentHashMap<>(this.ptms.size(),0.75F,cpu);



        long ri = 0L;
        for(Modification modification : this.ptms){


            if (modification.getModificationType() == ModificationType.modn_protein || modification.getModificationType() == ModificationType.modnaa_protein ||
                    modification.getModificationType() == ModificationType.modc_protein || modification.getModificationType() == ModificationType.modcaa_protein){
                continue;
            }

            if(modification.getMass() >= leftMass && modification.getMass() <= rightMass){
                ri = ri + 1L;
                ArrayList<ScoreResult> sList = new ArrayList<>();
                res.put(ri,sList);
                fixedThreadPool.execute(new ModPepQueryAndScoringWorker(ri,modification,spectrum, peptideIndexMap,res,targetPeptideSeq));
            }

        }

        fixedThreadPool.shutdown();

        try {
            fixedThreadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }


        return(res);

    }





    private boolean validation_status = false;
    public boolean get_validation_status(){
        return this.validation_status;
    }
    public synchronized void update_validation_status(boolean status){
        this.validation_status = status;
    }

    /**
     * For each target spectrum, find matched reference peptides with any modifications from UniMod.
     * @param spectrum_title Spectrum title
     * @param peptideIndexMap Peptide indexing hash map which contains all reference peptides. The key is peptide mass index (integer value). The value for each key is a PeptideIndex object.
     * @param targetPeptideSeq Target peptide sequence matched to the target spectrum.
     * @return A hash map which contains scoring result.
     */
    ConcurrentHashMap<Long, ArrayList<ScoreResult>> getCandidateRefPeptidesSingleThread(String spectrum_title,
                                                                                        HashMap<Long,ArrayList<JPeptide>> peptideIndexMap,
                                                                                        String targetPeptideSeq) {

        Spectrum spectrum = SpectraInput.spectraMap.get(spectrum_title);

        if(leftMass < getMinModificationMass()){
            leftMass = getMinModificationMass();
        }

        if(debug) {
            Cloger.getInstance().logger.info("left mass:" + leftMass);
            Cloger.getInstance().logger.info("left mass:" + rightMass);
        }
        double precursorMass = spectrum.getPrecursor().getMass(spectrum.getPrecursor().possibleCharges[0]);
        if(debug) {
            Cloger.getInstance().logger.info("precursorMass=" + precursorMass + ",leftMass=" + leftMass + ",rightMass=" + rightMass);
        }

        // HyperscoreMatch.generateFactorialValues(60);


        // peptide mass with modifications
        // The max and min delta mass (PTM mass) will be considered.
        long maxIntValue = Math.round(precursorMass-leftMass+2);
        long minIntValue = Math.round(precursorMass-rightMass-2);
        //System.out.println("left mass:" + leftMass);
        //System.out.println("left mass:" + rightMass);
        //System.out.println("precursorMass=" + precursorMass + ",leftMass=" + leftMass + ",rightMass=" + rightMass);

        int cpu = CParameter.cpu;

        if (cpu == 0) {
            cpu = Runtime.getRuntime().availableProcessors();
        }
        //System.out.println("CPU:" + cpu);
        ExecutorService fixedThreadPool = Executors.newFixedThreadPool(cpu);

        ConcurrentHashMap<Long,ArrayList<ScoreResult>> res = new ConcurrentHashMap<>(500,0.75F,cpu);

        //int ind = 1;

        for (long i = minIntValue; i <= maxIntValue; i++) {
            ArrayList<ScoreResult> sList = new ArrayList<>();
            res.put(i,sList);
            ModPepQueryAndScoringWorker modPepQueryAndScoringWorker = new ModPepQueryAndScoringWorker(i,spectrum, peptideIndexMap,res,targetPeptideSeq);
            boolean stop = modPepQueryAndScoringWorker.do_score();
            if(stop){
                break;
            }
        }

        return(res);

    }


    /**
    public String doPTMValidationMHC(String psmfile, ArrayList<PeptideInput> peptideInputs, String tagFile, ArrayList<String> fixMods, String outdir){

        long startTime = System.currentTimeMillis();

        // read PSM file and only perform PTM validation for confident PSMs
        BufferedReader bReader = null;
        try {
            bReader = new BufferedReader(new FileReader(new File(psmfile)));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        String line = null;
        try {
            line = bReader.readLine().trim();
        } catch (IOException e) {
            e.printStackTrace();
        }
        String headLine[] = line.split("\t");
        HashMap<String,Integer> headMap = new HashMap<>();
        for(int i=0;i<headLine.length;i++){
            headMap.put(headLine[i],i);
        }


        HashMap<String,Spectrum> spectraMap = new HashMap<>();

        HashMap<String,String> spectra2pep = new HashMap<>();

        try {
            while((line = bReader.readLine())!=null){
                line = line.trim();
                String d[] = line.split("\t");
                String spectrum_title  = d[headMap.get("spectrum_title")];
                String peptideSequence = d[headMap.get("peptide")];

                double pvalue = Double.parseDouble(d[headMap.get("pvalue")]);
                int rank = Integer.parseInt(d[headMap.get("rank")]);
                int n_db = Integer.parseInt(d[headMap.get("n_db")]);

                boolean psm_validation_passed = CParameter.passed_psm_validation(peptideSequence.length(),
                        pvalue,n_db,rank);

                if(psm_validation_passed){
                    spectraMap.put(spectrum_title,new Spectrum());
                    spectra2pep.put(spectrum_title,peptideSequence);
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

        String outfile = outdir + "/ptm.txt";

        if(spectraMap.size()>=1) {

            for (PeptideInput peptideInput : peptideInputs) {
                if (peptideInput.getOutLines().size() >= 1) {
                    //String peptideSequence = peptideInput.jPeptide.peptide.getSequence();
                    //for (MSnSpectrum spectrum : peptideInput.jPeptide.mSnSpectrums) {


                    for (JPeptide jPeptide : peptideInput.getPtmIsoforms()) {
                        //for (MSnSpectrum spectrum : jPeptide.mSnSpectrums) {
                        for (String title : jPeptide.spectraIndexs) {
                            Spectrum spectrum = SpectraInput.spectraMap.get(title);
                            if (spectraMap.containsKey(spectrum.getSpectrumTitle())) {
                                spectraMap.put(spectrum.getSpectrumTitle(), spectrum);
                            }
                        }
                    }
                }
            }

            Cloger.getInstance().logger.info("The number of spectra that need to be validated by PTM searching: " + spectraMap.size());

            //try {
            //    importPTMsFromUnimod();
            //} catch (IOException e) {
            //    e.printStackTrace();
            //}

            // HashMap<Long, PeptideIndex> peptideIndexMap = buildPeptideDB(db, fixMods);

            HashSet<String> spectraIDs = new HashSet<>(spectraMap.keySet());
            HashMap<String, HashMap<Long, PeptideIndex>> peptideIndexMap = new HashMap<>();
            try {
                peptideIndexMap = buildPeptideDB( tagFile, spectraIDs, 0 );
            } catch (IOException e) {
                e.printStackTrace();
            }


            BufferedWriter bWriter = null;
            try {
                bWriter = new BufferedWriter(new FileWriter(new File(outfile)));
            } catch (IOException e) {
                e.printStackTrace();
            }
            if(CParameter.scoreMethod==0) {
                try {
                    bWriter.write("spectrum_title\tpeptide\tcharge\texp_mass\tpep_mass\tmodification\tscore1\tscore2\n");
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }else{
                try {
                    bWriter.write("spectrum_title\tpeptide\tcharge\texp_mass\tpep_mass\tmodification\tscore\n");
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

            // add fixed modification into peptide
            ModPepQueryAndScoringWorker.fixMods = fixMods;

            // Perform the search for each target spectrum.
            for (String title : spectraMap.keySet()) {

                if(!peptideIndexMap.containsKey(title)){
                    continue;
                }

                ConcurrentHashMap<Long, ArrayList<ScoreResult>> res = this.getCandidateRefPeptides(spectraMap.get(title), peptideIndexMap.get(title), spectra2pep.get(title), CParameter.cpu);

                for(long l : res.keySet()){
                    for(ScoreResult scoreResult: res.get(l)){
                        StringBuilder stringBuilder = new StringBuilder();
                        Peptide peptide = scoreResult.peptide;
                        double exp_mass = spectraMap.get(title).getPrecursor().getMass(spectraMap.get(title).getPrecursor().possibleCharges[0]);
                        int charge = spectraMap.get(title).getPrecursor().possibleCharges[0];
                        double pep_mass = JPeptide.getMass(peptide);
                        if(CParameter.scoreMethod==0) {
                            stringBuilder.append(title).append("\t")
                                    .append(peptide.getSequence()).append("\t")
                                    .append(charge).append("\t")
                                    .append(getModificationString(peptide)).append("\t")
                                    .append(exp_mass).append("\t")
                                    .append(pep_mass).append("\t")
                                    .append(scoreResult.score).append("\t").append(scoreResult.score2).append("\n");
                        }else{
                            stringBuilder.append(title).append("\t")
                                    .append(peptide.getSequence()).append("\t")
                                    .append(charge).append("\t")
                                    .append(exp_mass).append("\t")
                                    .append(pep_mass).append("\t")
                                    .append(getModificationString(peptide)).append("\t").append(scoreResult.score).append("\n");
                        }

                        try {
                            bWriter.write(stringBuilder.toString());
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }
                }
            }

            try {
                bWriter.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        long endTime = System.currentTimeMillis();

        Cloger.getInstance().logger.info("Time used for PTM validation:" + (endTime - startTime) / 1000 + " s.");
        return(outfile);
    }**/



    public String doPTMValidation(String psmfile, ArrayList<PeptideInput> peptideInputs, String db, ArrayList<String> fixMods, String outdir){

        long startTime = System.currentTimeMillis();

        // read PSM file and only perform PTM validation for confident PSMs
        BufferedReader bReader = null;
        try {
            bReader = new BufferedReader(new FileReader(new File(psmfile)));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        String line = null;
        try {
            line = bReader.readLine().trim();
        } catch (IOException e) {
            e.printStackTrace();
        }
        String headLine[] = line.split("\t");
        HashMap<String,Integer> headMap = new HashMap<>();
        for(int i=0;i<headLine.length;i++){
            headMap.put(headLine[i],i);
        }


        HashMap<String,Spectrum> spectraMap = new HashMap<>();

        HashMap<String,String> spectra2pep = new HashMap<>();

        // store the lowest score (min score) for each spectrum.
        // for each spectrum, only need to save PSMs with scores >= min score in the unrestricted modification searching.
        // this can save space
        HashMap<String, Double> spectra2score = new HashMap<>();

        // save target peptides which are the peptides needed to be validated.
        HashSet<String> target_peptides = new HashSet<>();

        try {
            while((line = bReader.readLine())!=null){
                line = line.trim();
                String d[] = line.split("\t");
                String spectrum_title  = d[headMap.get("spectrum_title")];
                String peptideSequence = d[headMap.get("peptide")];
                target_peptides.add(peptideSequence.replaceAll("I","L"));

                double pvalue = Double.parseDouble(d[headMap.get("pvalue")]);
                int rank = Integer.parseInt(d[headMap.get("rank")]);
                int n_db = Integer.parseInt(d[headMap.get("n_db")]);
                double score = Double.parseDouble(d[headMap.get("score")]);

                boolean psm_validation_passed = CParameter.passed_psm_validation(peptideSequence.length(),
                        pvalue,n_db,rank);

                if(psm_validation_passed){
                    spectraMap.put(spectrum_title,new Spectrum());
                    spectra2pep.put(spectrum_title,peptideSequence);

                    // save the min score for each spectrum needed to validate.
                    if(spectra2score.containsKey(spectrum_title)){
                        if(spectra2score.get(spectrum_title) > score){
                            spectra2score.put(spectrum_title,score);
                        }
                    }else{
                        spectra2score.put(spectrum_title,score);
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

        String outfile = outdir + "/ptm.txt";

        if(spectraMap.size()>=1) {

            Cloger.getInstance().logger.info("The number of spectra that need to be validated by PTM searching: " + spectraMap.size());

            //try {
            //    ModificationDB.importPTMsFromUnimod();
            //} catch (IOException e) {
            //    e.printStackTrace();
            //}

            //HashMap<Long, PeptideIndex> peptideIndexMap = buildPeptideDB(db, fixMods, target_peptides);
            //HashMap<Long, PeptideIndex> peptideIndexMap = buildPeptideDB(db, target_peptides);
            //HashMap<Long, ArrayList<JPeptide>> peptideIndexMap = DatabaseInput.generate_reference_peptide_index(db,target_peptides,1.0);
            //HashMap<Long, ArrayList<JPeptide>> peptideIndexMap = DatabaseInput.generate_reference_peptide_index(db,target_peptides,10.0);
            if(ModificationDB.getInstance().get_peptide_index().size()<=0){
                ModificationDB.getInstance().add_peptide_index(DatabaseInput.generate_reference_peptide_index(db,target_peptides,10.0));
            }



            BufferedWriter bWriter = null;
            try {
                bWriter = new BufferedWriter(new FileWriter(new File(outfile)));
            } catch (IOException e) {
                e.printStackTrace();
            }
            if(CParameter.scoreMethod==0) {
                try {
                    bWriter.write("spectrum_title\tpeptide\tcharge\texp_mass\tpep_mass\ttol_ppm\ttol_da\tisotope_error\tmodification\tscore1\tscore2\n");
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }else{
                try {
                    bWriter.write("spectrum_title\tpeptide\tcharge\texp_mass\tpep_mass\ttol_ppm\ttol_da\tisotope_error\tmodification\tscore\n");
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

            // add fixed modification into peptide
            ModPepQueryAndScoringWorker.fixMods = fixMods;
            ModPepQueryAndScoringWorker.spectra2score = spectra2score;

            // Perform the search and scoring for each target spectrum.
            for (String title : spectraMap.keySet()) {
                this.update_validation_status(false);
                //ConcurrentHashMap<Long, ArrayList<ScoreResult>> res = this.getCandidateRefPeptides(spectraMap.get(title), peptideIndexMap, spectra2pep.get(title));
                ConcurrentHashMap<Long, ArrayList<ScoreResult>> res = this.getCandidateRefPeptidesFast(SpectraInput.spectraMap.get(title), ModificationDB.getInstance().get_peptide_index(), spectra2pep.get(title));
                this.save(bWriter, res, spectra2score, SpectraInput.spectraMap.get(title));

            }

            try {
                bWriter.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        long endTime = System.currentTimeMillis();
        ModificationDB.getInstance().get_peptide_index().clear();

        ModPepQueryAndScoringWorker.spectra2score = null;

        Cloger.getInstance().logger.info("Time used for PTM validation:" + (endTime - startTime) / 1000 + " s.");
        return(outfile);
    }


    public String doPTMValidationFast(String psmfile, ArrayList<PeptideInput> peptideInputs, String db, ArrayList<String> fixMods, String outdir){

        long startTime = System.currentTimeMillis();

        // read PSM file and only perform PTM validation for confident PSMs
        BufferedReader bReader = null;
        try {
            bReader = new BufferedReader(new FileReader(new File(psmfile)));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        String line = null;
        try {
            line = bReader.readLine().trim();
        } catch (IOException e) {
            e.printStackTrace();
        }
        String headLine[] = line.split("\t");
        HashMap<String,Integer> headMap = new HashMap<>();
        for(int i=0;i<headLine.length;i++){
            headMap.put(headLine[i],i);
        }


        HashMap<String,Spectrum> spectraMap = new HashMap<>();

        HashMap<String,String> spectra2pep = new HashMap<>();

        // store the lowest score (min score) for each spectrum.
        // for each spectrum, only need to save PSMs with scores >= min score in the unrestricted modification searching.
        // this can save space
        HashMap<String, Double> spectra2score = new HashMap<>();

        // save target peptides which are the peptides needed to be validated.
        HashSet<String> target_peptides = new HashSet<>();

        try {
            while((line = bReader.readLine())!=null){
                line = line.trim();
                String d[] = line.split("\t");
                String spectrum_title  = d[headMap.get("spectrum_title")];
                String peptideSequence = d[headMap.get("peptide")];
                target_peptides.add(peptideSequence.replaceAll("I","L"));

                double pvalue = Double.parseDouble(d[headMap.get("pvalue")]);
                int rank = Integer.parseInt(d[headMap.get("rank")]);
                int n_db = Integer.parseInt(d[headMap.get("n_db")]);
                double score = Double.parseDouble(d[headMap.get("score")]);

                boolean psm_validation_passed = CParameter.passed_psm_validation(peptideSequence.length(),
                        pvalue,n_db,rank);

                if(psm_validation_passed){
                    spectraMap.put(spectrum_title,new Spectrum());
                    spectra2pep.put(spectrum_title,peptideSequence);

                    // save the min score for each spectrum needed to validate.
                    if(spectra2score.containsKey(spectrum_title)){
                        //if(spectra2score.get(spectrum_title) > score){
                        if(spectra2score.get(spectrum_title) < score){
                            spectra2score.put(spectrum_title,score);
                        }
                    }else{
                        spectra2score.put(spectrum_title,score);
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

        String outfile = outdir + "/ptm.txt";

        if(spectraMap.size()>=1) {

            for (PeptideInput peptideInput : peptideInputs) {
                if (peptideInput.getOutLines().size() >= 1) {
                    //String peptideSequence = peptideInput.jPeptide.peptide.getSequence();
                    //for (MSnSpectrum spectrum : peptideInput.jPeptide.mSnSpectrums) {
                    /**
                     for (String title : peptideInput.jPeptide.spectraIndexs) {
                     Spectrum spectrum = SpectraInput.spectraMap.get(title);
                     if (spectraMap.containsKey(spectrum.getSpectrumTitle())) {
                     spectraMap.put(spectrum.getSpectrumTitle(), spectrum);
                     }
                     }**/

                    for (JPeptide jPeptide : peptideInput.getPtmIsoforms()) {
                        //for (MSnSpectrum spectrum : jPeptide.mSnSpectrums) {
                        for (String title : jPeptide.spectraIndexs) {
                            Spectrum spectrum = SpectraInput.spectraMap.get(title);
                            if (spectraMap.containsKey(spectrum.getSpectrumTitle())) {
                                spectraMap.put(spectrum.getSpectrumTitle(), spectrum);
                            }
                        }
                    }
                }
            }

            Cloger.getInstance().logger.info("The number of spectra that need to be validated by PTM searching: " + spectraMap.size());

            //try {
            //    ModificationDB.importPTMsFromUnimod();
            //} catch (IOException e) {
            //    e.printStackTrace();
            //}

            //HashMap<Long, PeptideIndex> peptideIndexMap = buildPeptideDB(db, fixMods, target_peptides);
            //HashMap<Long, PeptideIndex> peptideIndexMap = buildPeptideDB(db, target_peptides);
            HashMap<Long, ArrayList<JPeptide>> peptideIndexMap = DatabaseInput.generate_reference_peptide_index(db,target_peptides,1.0);



            BufferedWriter bWriter = null;
            try {
                bWriter = new BufferedWriter(new FileWriter(new File(outfile)));
            } catch (IOException e) {
                e.printStackTrace();
            }
            if(CParameter.scoreMethod==0) {
                try {
                    bWriter.write("spectrum_title\tpeptide\tcharge\texp_mass\tpep_mass\ttol_ppm\ttol_da\tisotope_error\tmodification\tscore1\tscore2\n");
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }else{
                try {
                    bWriter.write("spectrum_title\tpeptide\tcharge\texp_mass\tpep_mass\ttol_ppm\ttol_da\tisotope_error\tmodification\tscore\n");
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

            // add fixed modification into peptide
            ModPepQueryAndScoringWorker.fixMods = fixMods;


            int cpu = 1;
            if(CParameter.cpu <= 0){
                cpu = Runtime.getRuntime().availableProcessors();
            }else{
                cpu = CParameter.cpu;
            }
            ExecutorService fixedThreadPool = Executors.newFixedThreadPool(cpu);

            SpectrumMatchWorker.spectra2score = spectra2score;
            SpectrumMatchWorker.peptideIndexMap = peptideIndexMap;
            ModPepQueryAndScoringWorker.spectra2score = spectra2score;
            // Perform the search and scoring for each target spectrum.
            for (String title : spectraMap.keySet()) {
                fixedThreadPool.execute(new SpectrumMatchWorker(title,spectra2pep.get(title),bWriter));
            }

            fixedThreadPool.shutdown();

            try {
                fixedThreadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }


            try {
                bWriter.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        long endTime = System.currentTimeMillis();

        Cloger.getInstance().logger.info("Time used for PTM validation:" + (endTime - startTime) / 1000 + " s.");
        return(outfile);
    }

    public synchronized void save(BufferedWriter bWriter,ConcurrentHashMap<Long, ArrayList<ScoreResult>> res, HashMap<String, Double> spectra2score, Spectrum spectrum){
        if(res.size()>=1) {
            // This is used for saving the best matching from the reference database searching
            // for the target spectrum
            ScoreResult bestRefScoreResult = new ScoreResult();
            // This is used to store already saved PSMs from reference database searching in this step.
            // The purpose of this is to avoid duplicated store of best PSM from this reference database searching
            HashSet<String> savedRefPSMs = new HashSet<>();

            for (long l : res.keySet()) {
                for (ScoreResult scoreResult : res.get(l)) {

                    if (scoreResult.score >= spectra2score.get(spectrum.getSpectrumTitle())) {
                        String out_line = get_psm_output(scoreResult, spectrum);
                        savedRefPSMs.add(out_line);
                        try {
                            bWriter.write(out_line);
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }

                    // We need to export the best PSM the reference database searching for
                    // calculating mod_delta_score
                    if (scoreResult.score > bestRefScoreResult.score) {
                        bestRefScoreResult = scoreResult;
                    }
                }
            }

            if (bestRefScoreResult.score > 0) {
                // save the best matching from the reference database searching for the target spectrum
                String bestRefPSM = get_psm_output(bestRefScoreResult, spectrum);
                if (!savedRefPSMs.contains(bestRefPSM)) {
                    try {
                        bWriter.write(bestRefPSM);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }
        }
    }

    public String getModificationString(Peptide peptide){
        String mod = "-";
        ModificationMatch modificationMatchs[] = peptide.getVariableModifications();
        if(modificationMatchs !=null && modificationMatchs.length>=1){
            String d[] = new String[modificationMatchs.length];
            for(int i=0;i<d.length;i++){
                Modification ptm = ptmFactory.getModification(modificationMatchs[i].getModification());
                d[i] = ptm.getName()+"@"+modificationMatchs[i].getSite()+"["+String.format("%.4f",ptm.getMass())+"]";
            }
            mod = StringUtils.join(d,';');
        }

        return(mod);
    }


    private SequenceProvider sequenceProvider = new SingleProteinSequenceProvider();
    private SequenceMatchingParameters modificationsSequenceMatchingParameters = SequenceMatchingParameters.getDefaultSequenceMatching();

    /**
     * Add fixed modifications to the peptide.
     * The input peptide object will be changed.
     * @param fixMods Fixed modification
     * @param peptide A Peptide object
     */
    public void addFixedModification(ArrayList<String> fixMods, Peptide peptide){

        for (String mod : fixMods) {
            Modification ptm = ptmFactory.getModification(mod);
            //ArrayList<Integer> possibleSites = peptide.getPotentialModificationSitesNoCombination(ptm, peptide.getSequence(), 0);
            int possibleSites[] = ModificationUtils.getPossibleModificationSites(peptide, ptm, sequenceProvider, modificationsSequenceMatchingParameters);

            for (Integer k : possibleSites) {
                peptide.addVariableModification(new ModificationMatch(ptm.getName(),k));
            }
        }


    }


    /**
     * Generate output for a PSM
     * @param scoreResult A ScoreResult object which stores score and peptide information used in matching
     * @param spectrum Spectrum object
     * @return A string which contains information needed to be saved to output file
     */
    public String get_psm_output(ScoreResult scoreResult, Spectrum spectrum){
        StringBuilder stringBuilder = new StringBuilder();
        Peptide peptide = scoreResult.peptide;
        double exp_mass = spectrum.getPrecursor().getMass(spectrum.getPrecursor().possibleCharges[0]);
        int charge = spectrum.getPrecursor().possibleCharges[0];
        double pep_mass = JPeptide.getMass(peptide);
        ArrayList<Double> tol_res = CParameter.get_mass_error(exp_mass, pep_mass);
        if (CParameter.scoreMethod == 0) {
            stringBuilder.append(spectrum.getSpectrumTitle()).append("\t")
                    .append(peptide.getSequence()).append("\t")
                    .append(charge).append("\t")
                    .append(exp_mass).append("\t")
                    .append(pep_mass).append("\t")
                    .append(tol_res.get(0)).append("\t")
                    .append(tol_res.get(1)).append("\t")
                    .append(tol_res.get(2)).append("\t")
                    .append(getModificationString(peptide)).append("\t")
                    .append(scoreResult.score).append("\t").append(scoreResult.score2).append("\n");
        } else {
            stringBuilder.append(spectrum.getSpectrumTitle()).append("\t")
                    .append(peptide.getSequence()).append("\t")
                    .append(charge).append("\t")
                    .append(exp_mass).append("\t")
                    .append(pep_mass).append("\t")
                    .append(tol_res.get(0)).append("\t")
                    .append(tol_res.get(1)).append("\t")
                    .append(tol_res.get(2)).append("\t")
                    .append(getModificationString(peptide)).append("\t").append(scoreResult.score).append("\n");
        }
        return stringBuilder.toString();
    }

}
