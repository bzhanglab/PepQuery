package main.java.OpenModificationSearch;

import com.compomics.util.experiment.biology.modifications.Modification;
import com.compomics.util.experiment.biology.modifications.ModificationType;
import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import main.java.PSMMatch.JPeptide;
import main.java.pg.CParameter;
import main.java.PSMMatch.ScoreResult;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ConcurrentHashMap;

import static main.java.pg.PeptideSearchMT.scorePeptide2Spectrum;

public final class ModPepQueryAndScoringWorker implements Runnable{


    /**
     * The best score for each spectrum matched to targeted peptides
     */
    public static HashMap<String,Double> spectra2score = new HashMap<>();


    public static ArrayList<String> fixMods;
    // private int ind = 0;
    /**
     * peptide index in integer value format
     */
    private long intValue;

    private Modification modification;

    private Spectrum spectrum = null;
    private HashMap<Long, ArrayList<JPeptide>> peptideIndexMap = null;

    private ConcurrentHashMap<Long,ArrayList<ScoreResult>> res = null;
    private String targetPeptideSeq  = null;

    public ModPepQueryAndScoringWorker(long i, Spectrum spectrumObj, HashMap<Long, ArrayList<JPeptide>> peptideIndMap, ConcurrentHashMap<Long,ArrayList<ScoreResult>> resMap) {
        this.intValue = i;
        this.spectrum = spectrumObj;
        this.peptideIndexMap = peptideIndMap;
        this.res = resMap;
    }

    public ModPepQueryAndScoringWorker(long i, Spectrum spectrumObj, HashMap<Long, ArrayList<JPeptide>> peptideIndMap, ConcurrentHashMap<Long,ArrayList<ScoreResult>> resMap, String targetPepSeq) {
        this.intValue = i;
        this.spectrum = spectrumObj;
        this.peptideIndexMap = peptideIndMap;
        this.res = resMap;
        this.targetPeptideSeq = targetPepSeq;
    }

    public ModPepQueryAndScoringWorker(long i, Modification modification, Spectrum spectrum, HashMap<Long, ArrayList<JPeptide>> peptideIndMap, ConcurrentHashMap<Long,ArrayList<ScoreResult>> resMap, String targetPepSeq) {
        this.intValue = i;
        this.modification = modification;
        this.spectrum = spectrum;
        this.peptideIndexMap = peptideIndMap;
        this.res = resMap;
        this.targetPeptideSeq = targetPepSeq;
    }

    @Override
    public void run() {
        if(!ModificationDB.getInstance().get_validation_status()){
            if(this.modification!=null){
                this.do_score(modification);
            }else {
                this.do_score();
            }
        }
    }

    public boolean do_score(){
        // save the best match from modification searching
        // double bestPSMscore = 0.0;
        // Peptide bestPeptide = new Peptide();
        // Only consider peptides belong to peptide index - intValue
        boolean stop = false;
        if(peptideIndexMap.containsKey(intValue)){
            ArrayList<JPeptide> peptides = peptideIndexMap.get(intValue);
            for(int k=0;k<peptides.size();k++){

                if(ModificationDB.getInstance().get_validation_status()){
                    return true;
                }

                /**
                if(this.targetPeptideSeq != null){
                    // calculate the similarity of the query peptide and target peptide.
                    if(!calcSimilarityOfPeptides(this.targetPeptideSeq, peptides.get(k).peptide.getSequence())){
                        continue;
                    }
                }**/

                ArrayList<Peptide> modPeptides = ModificationDB.getInstance().getCandidatePTMs(spectrum,peptides.get(k).peptide);
                for(Peptide pep: modPeptides){
                    if(ModificationDB.getInstance().get_validation_status()){
                        return true;
                    }
                    ArrayList<Double> tol_res = CParameter.get_mass_error(spectrum,pep);
                    // we need to check this condition since the peptides got from getCandidatePTMs may not fit this condition.
                    if(tol_res.size()>0) {
                        ScoreResult scoreResult = scorePeptide2Spectrum(pep, spectrum);

                        // save the best match from modification searching
                        // if(bestPSMscore < scoreResult.score){
                        //    bestPSMscore = scoreResult.score;
                        //    bestPeptide = pep;
                        //}
                        if (scoreResult.score >= CParameter.minScore) {
                            // Only output PSMs with score >= minScore
                            scoreResult.peptide = pep;
                            this.res.get(intValue).add(scoreResult);
                        }

                        if(CParameter.fast_model){
                            if(spectra2score.containsKey(spectrum.getSpectrumTitle())){
                                if(CParameter.unrestrictedSearchWithEqualAndBetterScore) {
                                    if (scoreResult.score >= spectra2score.get(spectrum.getSpectrumTitle())) {
                                        stop = true;
                                        break;
                                    }
                                }else{
                                    if (scoreResult.score > spectra2score.get(spectrum.getSpectrumTitle())) {
                                        stop = true;
                                        break;
                                    }
                                }
                            }
                        }
                        //System.out.println(Thread.currentThread().getName()+"\t"+ind+"\t"+pep.getSequenceWithLowerCasePtms()+"\t"+pep.getMass()+"\t"+spectrum.getPrecursor().getMass(spectrum.getPrecursor().getPossibleCharges().get(0).value)+"\t"+scoreResult.score);
                    }

                }
                if(stop){
                    ModificationDB.getInstance().update_validation_status(true);
                    break;
                }
            }
        }
        return stop;
    }


    public boolean do_score(Modification modification){

        String mod_aa = get_aa_of_modification(modification);
        // save the best match from modification searching
        // double bestPSMscore = 0.0;
        // Peptide bestPeptide = new Peptide();
        // Only consider peptides belong to peptide index - intValue
        boolean better_match_found = false;
        // This is the actual mass of spectrum from MS/MS data
        double precursorMass = spectrum.getPrecursor().getMass(spectrum.getPrecursor().possibleCharges[0]);
        // isotope_masses store mono mass for precursor
        ArrayList<Double> isotope_masses = CParameter.get_precursor_mass_with_isotope_error(precursorMass);
        for(double precursor_mono_mass: isotope_masses) {
            if (ModificationDB.getInstance().get_validation_status()) {
                return true;
            }
            double delta_mass = precursor_mono_mass - modification.getMass();
            long peptideInd = Math.round(10.0 * delta_mass);
            // peptide index must be already sorted.
            boolean spectra_matched = false;
            boolean exceed_tol_limit = false;
            for(long peptide_index=peptideInd-1;peptide_index<=peptideInd+1;peptide_index++) {
                if (ModificationDB.getInstance().get_validation_status()) {
                    return true;
                }
                if (peptideIndexMap.containsKey(peptide_index)) {
                    ArrayList<JPeptide> peptides = peptideIndexMap.get(peptide_index);

                    for (int k = 0; k < peptides.size(); k++) {

                        if (ModificationDB.getInstance().get_validation_status()) {
                            return true;
                        }

                        double pep_mass = peptides.get(k).getMass() + modification.getMass();
                        // ArrayList<Double> tol_res = CParameter.get_mass_error(precursor_mono_mass, pep_mass);
                        // we need to check this condition since the peptides got from getCandidatePTMs may not fit this condition.
                        // if (tol_res.size() <= 0) {
                        //     continue;
                        //}

                        double del = Math.abs(CParameter.get_mass_error(precursor_mono_mass, pep_mass, CParameter.tolu));

                        if (del <= CParameter.tol) {
                            spectra_matched = true;
                        }else{
                            if(spectra_matched){
                                // Since the peptide index is already sorted based on peptide mass from small to large,
                                // if there is already spectra matched and the current match has delta mass (del) > CParameter.tol,
                                // we no longer need to compare for the remaining peptides.
                                exceed_tol_limit = true;
                                break;
                            }
                            continue;
                        }

                        // if the modification occurs on an amino acid and the peptide sequence doesn't contain this
                        // amino acid, will ignore this peptide.
                        if(mod_aa!=null && !peptides.get(k).peptide.getSequence().contains(mod_aa)){
                            continue;
                        }


                        /**
                        if (this.targetPeptideSeq != null) {
                            // calculate the similarity of the query peptide and target peptide.
                            if (!calcSimilarityOfPeptides(this.targetPeptideSeq, peptides.get(k).peptide.getSequence())) {
                                continue;
                            }
                        }**/

                        ArrayList<Peptide> modPeptides = ModificationDB.getInstance().getCandidatePTMs(peptides.get(k).peptide, modification);
                        for (Peptide pep : modPeptides) {
                            if (ModificationDB.getInstance().get_validation_status()) {
                                return true;
                            }

                            ScoreResult scoreResult = scorePeptide2Spectrum(pep, spectrum);

                            // save the best match from modification searching
                            // if(bestPSMscore < scoreResult.score){
                            //    bestPSMscore = scoreResult.score;
                            //    bestPeptide = pep;
                            //}
                            if (scoreResult.score >= CParameter.minScore) {
                                // Only output PSMs with score >= minScore
                                scoreResult.peptide = pep;
                                this.res.get(intValue).add(scoreResult);
                            }

                            if (CParameter.fast_model) {
                                if (spectra2score.containsKey(spectrum.getSpectrumTitle())) {
                                    if (CParameter.unrestrictedSearchWithEqualAndBetterScore) {
                                        if (scoreResult.score >= spectra2score.get(spectrum.getSpectrumTitle())) {
                                            better_match_found = true;
                                            break;
                                        }
                                    } else {
                                        if (scoreResult.score > spectra2score.get(spectrum.getSpectrumTitle())) {
                                            better_match_found = true;
                                            break;
                                        }
                                    }
                                }
                            }

                        }
                        if (better_match_found) {
                            ModificationDB.getInstance().update_validation_status(true);
                            break;
                        }

                    }
                    if (better_match_found) {
                        ModificationDB.getInstance().update_validation_status(true);
                        break;
                    }
                    if(exceed_tol_limit){
                        // This just means we don't need to compare the remaining peptides due to tol limit issue.
                        // It doesn't mean there is a better match found, so we cannot use:
                        // ModificationDB.getInstance().update_validation_status(true);
                        break;
                    }
                }
            }
        }
        return better_match_found;
    }

    /**
    // 2022-04-01
    public void run_bak() {

        // save the best match from modification searching
        // double bestPSMscore = 0.0;
        // Peptide bestPeptide = new Peptide();
        // Only consider peptides belong to peptide index - intValue
        if(peptideIndexMap.containsKey(intValue)){
            ArrayList<String> peptideSequences = peptideIndexMap.get(intValue).peptideSequences;
            ArrayList<Double> peptideMasses = peptideIndexMap.get(intValue).peptideMasses;
            for(int k=0;k<peptideMasses.size();k++){

                if(this.targetPeptideSeq != null){
                    // calculate the similarity of the query peptide and target peptide.
                    if(!calcSimilarityOfPeptides(this.targetPeptideSeq, peptideSequences.get(k))){
                        continue;
                    }
                }

                Peptide peptide = new Peptide(peptideSequences.get(k));
                // must consider fixed modification first, otherwise it's difficult to find potential variable modification
                // candidate
                ModificationDB.getInstance().addFixedModification(fixMods,peptide);
                ArrayList<Peptide> modPeptides = ModificationDB.getInstance().getCandidatePTMs(spectrum,peptide);
                for(Peptide pep: modPeptides){
                    ArrayList<Double> tol_res = CParameter.get_mass_error(spectrum,pep);
                    // we need to check this condition since the peptides got from getCandidatePTMs may not fit this condition.
                    if(tol_res.size()>0) {
                        ScoreResult scoreResult = scorePeptide2Spectrum(pep, spectrum);

                        // save the best match from modification searching
                        // if(bestPSMscore < scoreResult.score){
                        //    bestPSMscore = scoreResult.score;
                        //    bestPeptide = pep;
                        //}
                        if (scoreResult.score >= CParameter.minScore) {
                            // Only output PSMs with score >= minScore
                            scoreResult.peptide = pep;
                            this.res.get(intValue).add(scoreResult);
                        }
                        //System.out.println(Thread.currentThread().getName()+"\t"+ind+"\t"+pep.getSequenceWithLowerCasePtms()+"\t"+pep.getMass()+"\t"+spectrum.getPrecursor().getMass(spectrum.getPrecursor().getPossibleCharges().get(0).value)+"\t"+scoreResult.score);
                    }

                }
            }
        }

    }**/


    public boolean calcSimilarityOfPeptides(String seq1, String seq2){
        String d1[] = seq1.split("");
        String d2[] = seq2.split("");
        HashSet<String> h1 = new HashSet<>();

        boolean yes = false;
        for(int i=0;i<(d1.length-2);i++){
            h1.add(d1[i]+d1[i+1]);
        }

        for(int i=0;i<(d2.length-2);i++){
            if(h1.contains(d2[i]+d2[i+1])){
                yes = true;
                break;
            }
        }

        return(yes);
    }


    /**
     * Only need to consider the case when the modification occurs on amino acid.
     * @param modification A Modification object.
     * @return null or amino acid string
      */
    public String get_aa_of_modification(Modification modification){
        if(modification.getModificationType() == ModificationType.modaa){
            ArrayList<Character> aa = modification.getPattern().getAminoAcidsAtTarget();
            if(aa.size()==1){
                return String.valueOf(aa.get(0));
            }else{
                //System.err.println("Modification error: multiple sites, "+modification.getName());
                return null;
            }
        }else{
            //System.err.println("Modification error: not aa modification, "+modification.getName());
            return null;
        }
    }


}


