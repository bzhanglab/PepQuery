package main.java.pg;


import com.compomics.util.experiment.biology.modifications.Modification;
import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.experiment.identification.protein_sequences.SingleProteinSequenceProvider;
import com.compomics.util.experiment.identification.utils.ModificationUtils;
import com.compomics.util.experiment.io.biology.protein.SequenceProvider;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import main.java.PSMMatch.JPeptide;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.paukov.combinatorics3.Generator;

import java.util.*;
import java.util.stream.Collectors;


public final class PeptideInput {

    //public JPeptide jPeptide = null;

    /**
     * Peptide sequence
     */
    public String peptideSequence = "";

    /**
     * All possible peptide forms for a given peptide sequence.
     */
    private ArrayList<JPeptide> isoforms = new ArrayList<>();

    /**
     * Saved all matched spectra to all the possible peptide forms of the target sequence.
     */
    private HashMap<String, Spectrum> spectra = new HashMap<>();

    public boolean hasSpectrum(String title){
        if(this.spectra.containsKey(title)){
            return true;
        }else{
            return false;
        }
    }

    public boolean hasSpectrum(Spectrum spectrum){
        if(this.spectra.containsKey(spectrum.getSpectrumTitle())){
            return true;
        }else{
            return false;
        }
    }

    public void addSpectrum(Spectrum spectrum){
        if(!this.spectra.containsKey(spectrum.getSpectrumTitle())){
            this.spectra.put(spectrum.getSpectrumTitle(),spectrum);
        }
    }

    public HashMap<String, Spectrum> getSpectra(){
        return this.spectra;
    }

    public void clearSpectra(){
        this.spectra.clear();
    }

    /**
     * Random peptides
     */
    public ArrayList<String> randomPeptides = new ArrayList<>();



    private static SequenceProvider sequenceProvider = new SingleProteinSequenceProvider();
    private static SequenceMatchingParameters modificationsSequenceMatchingParameters = SequenceMatchingParameters.getDefaultSequenceMatching();




    private ArrayList<String> outLines = new ArrayList<>();
    private ArrayList<String> outAdditionalLines = new ArrayList<>();
    private ArrayList<String> outPeakAnnotations = new ArrayList<>();

    public PeptideInput(String pep_seq){
        this.peptideSequence = pep_seq;
    }

    /**
     * Get all possible peptide forms for a given peptide sequence
     * @return ArrayList<JPeptide>
     */
    public ArrayList<JPeptide> getPtmIsoforms(){
        return(this.isoforms);
    }

    public void setPtmIsoforms(ArrayList<Peptide> ps){
        for(Peptide p:ps){
            JPeptide jPeptide = new JPeptide(p);
            this.isoforms.add(jPeptide);
        }
    }


    /**
     * For a given peptide sequence, generate all possible peptide forms considering defined fixed modification and
     * variable modification, as well as the max allowed number of variable modifications.
     * @param peptideSequence A peptide sequences
     * @return An ArrayList<Peptide> object.
     */
    public static ArrayList<Peptide> calcPeptideIsoforms(String peptideSequence){

        ArrayList<Modification> varMods = ModificationGear.getInstance().getPTMs(CParameter.varMods);
        int maxVarMods = CParameter.maxVarMods;

        Peptide peptide = new Peptide(peptideSequence);

        ArrayList<Peptide> peptides = new ArrayList<>();
        if(varMods.size()>=1) {
            ArrayList<JPTM> all_mod_sites = new ArrayList<>();
            // get all possible modifications sites for all variable modifications
            for (Modification ptm : varMods) {
                //ArrayList<Integer> possibleSites = pepObj.getPotentialModificationSitesNoCombination(ptm, pepObj.getSequence(), 0);
                int pSites[] = ModificationUtils.getPossibleModificationSites(peptide, ptm, sequenceProvider, modificationsSequenceMatchingParameters);

                if(pSites.length >= 1){
                    //ArrayList<Integer> possibleSites = new ArrayList<>();
                    for(int i : pSites){
                        all_mod_sites.add(new JPTM(i,ptm));
                    }
                }


            }
            if(all_mod_sites.size()>=1) {
                // >=1 possible modification site
                // the max number of modifications to consider
                int maxNumMods = Math.min(all_mod_sites.size(), maxVarMods);
                for(int k=1;k<=maxNumMods;k++){
                    for (List<JPTM> iCombination : Generator.combination(all_mod_sites).simple(k)) {
                        Peptide modPeptide = new Peptide(peptideSequence);
                        // same site only allows one variable or fixed modification. This is controlled by
                        // CParameter.maxModsPerAA.
                        // The max number of modifications occurring on any individual position
                        int maxModPerAA = 0;
                        // This is used to count how many modifications occurring on each position.
                        HashMap<Integer, Integer> pos2mods = new HashMap<>();
                        for (JPTM mod : iCombination) {
                            modPeptide.addVariableModification(new ModificationMatch(mod.modification.getName(), mod.pos));

                            if (pos2mods.containsKey(mod.pos)) {
                                pos2mods.put(mod.pos, pos2mods.get(mod.pos) + 1);
                            } else {
                                pos2mods.put(mod.pos, 1);
                            }
                            if (pos2mods.get(mod.pos) > maxModPerAA) {
                                maxModPerAA = pos2mods.get(mod.pos);
                            }
                        }

                        if (maxModPerAA <= CParameter.maxModsPerAA) {
                            // add fixed modifications
                            addFixedModification(modPeptide);
                            peptides.add(modPeptide);
                        }


                    }
                }
            }
        }

        // add peptide which only has fixed modifications
        Peptide peptide_no_mod = new Peptide(peptideSequence);
        addFixedModification(peptide_no_mod);
        peptides.add(peptide_no_mod);
        return(peptides);
    }


    /**
     * Add fixed modification into the peptide
     * @param modification Fixed modification index, the format is like : 1,2,3
     */
    public static void addFixedModification(String modification, Peptide peptide){
        // 0 means no fixed modification
        if(!modification.equalsIgnoreCase("0")) {
            String[] mods = modification.split(",");
            if (mods.length >= 1) {
                for (String mod : mods) {
                    Modification ptm = ModificationGear.getInstance().getPTM(Integer.valueOf(mod));
                    //ArrayList<Integer> possibleSites = this.jPeptide.peptide.getPotentialModificationSitesNoCombination(ptm, this.jPeptide.peptide.getSequence(), 0);

                    int possibleSites[] = ModificationUtils.getPossibleModificationSites(peptide, ptm, sequenceProvider, modificationsSequenceMatchingParameters);
                    // The location in the sequence. N-term modification is at index 0, C-term at index sequence length + 1, and other
                    // modifications at amino acid index starting from 1.
                    for (Integer k : possibleSites) {
                        peptide.addVariableModification(new ModificationMatch(ptm.getName(), k));
                    }
                }
            }
        }
    }

    /**
     * Add fixed modifications to a peptide.
     * Fixed modifications are added after adding variable modifications
     * @param peptide A Peptide object.
     */
    public static void addFixedModification(Peptide peptide){
        // 0 means no fixed modification
        if(!CParameter.fixMods.equalsIgnoreCase("0")) {
            ArrayList<Modification> fixedMods = ModificationGear.getInstance().getPTMs(CParameter.fixMods);
            if(fixedMods.size()>=1) {
                HashSet<Integer> varModSites = new HashSet<>();
                ModificationMatch []varModificationMatchs =  peptide.getVariableModifications();
                for(ModificationMatch match: varModificationMatchs){
                    varModSites.add(match.getSite());
                }
                for (Modification mod : fixedMods) {
                    int possibleSites[] = ModificationUtils.getPossibleModificationSites(peptide, mod, sequenceProvider, modificationsSequenceMatchingParameters);
                    // The location in the sequence. N-term modification is at index 0, C-term at index sequence length + 1, and other
                    // modifications at amino acid index starting from 1.
                    for (Integer k : possibleSites) {
                        if(!varModSites.contains(k)){
                            peptide.addVariableModification(new ModificationMatch(mod.getName(), k));
                        }

                    }
                }
            }
        }

    }

    public static void main(String[] args) {


    }


    /**
     * Generate random peptides for a target peptide
     * @param pepObj The object of target peptide
     * @param num The number of random peptides
     * @param fixCtermAA Fix the c-terminal amino acid
     * @return Random peptide sequences
     */
    public static ArrayList<String> generateRandomPeptides(Peptide pepObj, int num, boolean fixCtermAA){
        ArrayList<String> randomPeptides = new ArrayList<>();
        String peptideSequence = pepObj.getSequence();
        int peptideLength = peptideSequence.length();

        // No chage for c-term amino acid
        if(fixCtermAA){
            peptideLength = peptideLength - 1;
            // remove c-term amino acid
            peptideSequence = peptideSequence.substring(0,peptideLength);
        }

        // Calculate the max number of combinations
        HashMap<String,Integer> aaCount = new HashMap<>();
        for(int i=0;i<peptideLength;i++){
            String aa = String.valueOf(peptideSequence.charAt(i));
            if(aaCount.containsKey(aa)){
                aaCount.put(aa,1+aaCount.get(aa));
            }else{
                aaCount.put(aa,1);
            }
        }


        // sort aaCount, key=aa, value=int
        // Sort by key: low to high
        List<HashMap.Entry<String,Integer>> aaList=new ArrayList<>(aaCount.size());
        aaList.addAll(aaCount.entrySet());
        HashMapValueCompare vc = new HashMapValueCompare();
        aaList.sort(vc);

        // The total number of random peptides
        long nComb = 1;
        int rn = peptideLength;

        for(HashMap.Entry<String,Integer> it : aaList) {
            String aa = it.getKey();
            //
            int n = aaCount.get(aa);
            long c = CombinatoricsUtils.binomialCoefficient(rn, n);

            // if the number > Long.MAX_VALUE, the result may be a negative value (will cause problem)
            // so set a max value
            if(nComb >= 10000000){
                nComb = 10000000;
            }else{
                nComb = nComb * c;
            }
            rn = rn - n;
            //System.out.println(aa+"\t"+n+"\t"+c+"\t"+rn);
        }


        // generate random peptides
        HashSet<String> pSet = new HashSet<>(num);
        pSet.add(pepObj.getSequence());

        num = num>nComb? (int) nComb :num;
        //System.out.println("Combination: "+nComb+","+num);

        // Save the result for combinations from the first amino acid
        // lowest frequency.
        Set<Integer> pos1aa = new HashSet<>(peptideLength);
        for(int j=1;j<=peptideLength;j++){
            pos1aa.add(j);
        }
        List<List<Integer>> pp1aa = Generator.combination(pos1aa).simple(aaList.get(0).getValue()).stream().collect(Collectors.toList());

        for(int i=1;i<=num;i++){

            Set<Integer> pos = new HashSet<>(peptideLength);
            for(int j=1;j<=peptideLength;j++){
                pos.add(j);
            }
            String pep[] = new String[pepObj.getSequence().length()];
            for(int k=0;k<aaList.size();k++) {

                String aa = aaList.get(k).getKey();
                int  n = aaCount.get(aa);
                long c = CombinatoricsUtils.binomialCoefficient(pos.size(),n);
                int  r = new Random().nextInt((int)c);

                List<Integer> pp;
                if(k==0){
                    // speed up
                    pp = pp1aa.get(r);
                }else{
                    pp = Generator.combination(pos).simple(aaCount.get(aa)).stream().collect(Collectors.toList()).get(r);
                }

                for(int indpp : pp){
                    pep[indpp] = aa;
                }
                //pp.forEach(k -> pep[k]=aa);

                pos.removeAll(pp);
                //System.out.println(aa+"\t"+n+"\t"+c+"\t"+r+ "\t"+pp.stream().map(k -> String.valueOf(k)).collect(Collectors.joining("_")));

            }
            String rpep = StringUtils.join(pep,"")+pepObj.getSequence().charAt(peptideLength);

            if(pSet.contains(rpep)){
                continue;
            }else{
                pSet.add(rpep);
                randomPeptides.add(rpep);
            }

            //System.out.println(rpep);

        }

        return(randomPeptides);
    }



    /**
     * Generate random peptides for a target peptide. This is the version used currently.
     * @param pep_seq The object of target peptide
     * @param num The number of random peptides
     * @param fixCtermAA Fix the c-terminal amino acid
     * @return peptide sequences
     */
    public static ArrayList<String> generateRandomPeptidesFast(String pep_seq, int num, boolean fixCtermAA){
        ArrayList<String> randomPeptides = new ArrayList<>();
        String peptideSequence = pep_seq;
        int peptideLength = peptideSequence.length();

        // No chage for c-term amino acid
        if(fixCtermAA){
            peptideLength = peptideLength - 1;
            peptideSequence = peptideSequence.substring(0,peptideLength);
        }

        HashMap<String,Integer> aaCount = new HashMap<>();
        for(int i=0;i<peptideLength;i++){
            String aa = String.valueOf(peptideSequence.charAt(i));
            if(aaCount.containsKey(aa)){
                aaCount.put(aa,1+aaCount.get(aa));
            }else{
                aaCount.put(aa,1);
            }
        }

        List<HashMap.Entry<String,Integer>> aaList=new ArrayList<>(aaCount.size());
        aaList.addAll(aaCount.entrySet());
        HashMapValueCompare vc = new HashMapValueCompare();
        aaList.sort(vc);

        long nComb = 1;
        int rn = peptideLength;

        for(HashMap.Entry<String,Integer> it : aaList) {
            String aa = it.getKey();
            int n = aaCount.get(aa);
            long c = CombinatoricsUtils.binomialCoefficient(rn, n);

            if(nComb >= 10000000){
                nComb = 10000000;
            }else{
                nComb = nComb * c;
            }
            rn = rn - n;
            //System.out.println(aa+"\t"+n+"\t"+c+"\t"+rn);
        }

        HashSet<String> pSet = new HashSet<>(num);
        pSet.add(pep_seq);

        num = num>nComb? (int) nComb :num;
        //System.out.println("Combination: "+nComb+","+num);

        Set<Integer> pos1aa = new HashSet<>(peptideLength);
        for(int j=1;j<=peptideLength;j++){
            pos1aa.add(j);
        }
        List<List<Integer>> pp1aa = Generator.combination(pos1aa).simple(aaList.get(0).getValue()).stream().collect(Collectors.toList());

        // control the random with random seed
        Random random = new Random();
        random.setSeed(12457);

        for(int i=1;i<=num;i++){

            Set<Integer> pos = new HashSet<>(peptideLength);
            for(int j=1;j<=peptideLength;j++){
                pos.add(j);
            }
            String pep[] = new String[pep_seq.length()];
            for(int k=0;k<aaList.size();k++) {

                String aa = aaList.get(k).getKey();
                int  n = aaCount.get(aa);

                for(int rindex=0;rindex<n;rindex++){
                    long c = CombinatoricsUtils.binomialCoefficient(pos.size(),1);
                    int  r = random.nextInt((int)c);

                    int pp = Generator.combination(pos).simple(1).stream().collect(Collectors.toList()).get(r).get(0);

                    pep[pp] = aa;
                    //pp.forEach(k -> pep[k]=aa);

                    pos.remove(pp);
                }

                //System.out.println(aa+"\t"+n+"\t"+c+"\t"+r+ "\t"+pp.stream().map(k -> String.valueOf(k)).collect(Collectors.joining("_")));

            }
            String rpep = StringUtils.join(pep,"")+pep_seq.charAt(peptideLength);

            if(pSet.contains(rpep)){
                continue;
            }else{
                pSet.add(rpep);
                randomPeptides.add(rpep);
            }

            //System.out.println(rpep);

        }

        return(randomPeptides);
    }



    public void setRandomPeptides(ArrayList<String> rPeptides){
        this.randomPeptides = rPeptides;
    }


    public static Peptide getModificationPeptide(Peptide pepObj, String peptideSequence, String fixMods){

        // peptide object without any modification
        Peptide pep = new Peptide(peptideSequence);

        // already modified site information
        HashSet<Integer> usedPos = new HashSet<>();

        HashSet<Integer> fixedModSites = new HashSet<>();

        // add fixed modification
        HashSet<String> fixedModsMap = new HashSet<>();

        // this is used to control the variable site selection.
        Random random = new Random();
        random.setSeed(2022);

        // variable modification
        int n = pepObj.getVariableModifications().length;
        ModificationFactory ptmFactory = ModificationFactory.getInstance();
        if (n >= 1) {
            // modification data from base peptide object
            ModificationMatch modMatches[] = pepObj.getVariableModifications();
            for (ModificationMatch modificationMatch : modMatches) {

                // variable modification. The type of modification here is important.
                if (!fixedModsMap.contains(modificationMatch.getModification())) {
                    Modification ptm = ptmFactory.getInstance().getModification(modificationMatch.getModification());
                    // get the possible modification sites.
                    //ArrayList<Integer> possibleSites = pep.getPotentialModificationSitesNoCombination(ptm, pep.getSequence(), 0);
                    int pSites[] = ModificationUtils.getPossibleModificationSites(pep, ptm, sequenceProvider, modificationsSequenceMatchingParameters);

                    if(pSites.length >= 1){
                        ArrayList<Integer> possibleSites = new ArrayList<>();
                        for(int i: pSites){
                            possibleSites.add(i);
                        }
                        possibleSites.removeAll(fixedModSites);

                        // If the site has been considered
                        ArrayList<Integer> removeSites = new ArrayList<>();
                        for(int i : possibleSites){
                            if(usedPos.contains(i)){
                                removeSites.add(i);
                            }
                        }
                        if(removeSites.size()>=1){
                            possibleSites.removeAll(removeSites);
                        }

                        if(possibleSites.size()>=1) {

                            int index = random.nextInt(possibleSites.size());
                            // int index = ThreadLocalRandom.current().nextInt(possibleSites.size());
                            int ptmIndex = possibleSites.get(index);
                            // ptmIndex: the position of the modification in the sequence, 1 is the first residue
                            pep.addVariableModification(new ModificationMatch(ptm.getName(), ptmIndex));

                            usedPos.add(ptmIndex);
                        }
                    }



                }

            }
        }

        PeptideInput.addFixedModification(pep);
        return(pep);

    }

    public synchronized void addOutputLine(String line){
        this.outLines.add(line);
    }

    public ArrayList<String> getOutLines(){
        return(this.outLines);
    }

    public synchronized void addOutAdditionalLines(String line){
        this.outAdditionalLines.add(line);
    }

    public ArrayList<String> getOutAdditionalLines(){
        return(this.outAdditionalLines);
    }


    public synchronized void addOutPeakAnnotations(String line){
        this.outPeakAnnotations.add(line);
    }

    public ArrayList<String> getOutPeakAnnotations(){
        return(this.outPeakAnnotations);
    }




}
