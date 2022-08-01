package main.java.pg;

import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import main.java.PSMMatch.JPeptide;
import main.java.PSMMatch.ScoreResult;
import main.java.plot.PeakAnnotation;

import java.util.ArrayList;
import java.util.HashMap;

import static main.java.pg.PeptideSearchMT.scorePeptide2Spectrum;

public class DBSearchWorker implements  Runnable{


    public static HashMap<Long, ArrayList<JPeptide>> peptideIndexMap = new HashMap<>();

    private RefPepMatchResult refPepMatchResult;

    public DBSearchWorker(RefPepMatchResult refPepMatchResult){
        this.refPepMatchResult = refPepMatchResult;
    }

    @Override
    public void run() {
        Spectrum spectrum = SpectraInput.spectraMap.get(this.refPepMatchResult.spectrum_title);
        int charge = spectrum.getPrecursor().possibleCharges[0];
        double exp_mass = spectrum.getPrecursor().getMass(charge);

        ArrayList<Double> isotope_masses = CParameter.get_precursor_mass_with_isotope_error(exp_mass);
        // The number of matched peptide forms based on mass tol only.
        int n_matched_peptides = 0;
        int n_better = 0;
        ScoreResult bestRefScoreResult = new ScoreResult();
        ArrayList<String> out_data = new ArrayList<>(10);
        ArrayList<String> psm_annotation_data = new ArrayList<>(1);

        for(double mass:isotope_masses) {
            long va = Math.round(mass * 10.0);
            double peptideMass;
            double del;
            for (long i = (va - 1); i <= (va + 1); i++) {
                if (peptideIndexMap.containsKey(i)) {
                    // peptide index must be already sorted.
                    boolean spectra_matched = false;
                    for (JPeptide pep : peptideIndexMap.get(i)) {
                        peptideMass = pep.getMass();
                        del = Math.abs(peptideMass - mass);
                        if (CParameter.tolu.equalsIgnoreCase("ppm")) {
                            del = (1.0e6) * 1.0 * del / peptideMass;
                        }
                        if (del <= CParameter.tol) {
                            spectra_matched = true;
                            n_matched_peptides++;

                            // scoring
                            ScoreResult tmpScoreResult = scorePeptide2Spectrum(pep.peptide, spectrum);
                            //tmpScoreResult.peptide = ePeptide.peptide;
                            tmpScoreResult.jPeptide = pep;
                            tmpScoreResult.spectrum_title = spectrum.getSpectrumTitle();
                            double tmpScore = tmpScoreResult.score;

                            // save the best matching from the reference database searching for the target spectrum
                            if(bestRefScoreResult.score < tmpScore){
                                bestRefScoreResult = tmpScoreResult;
                            }

                            if (this.refPepMatchResult.score <= tmpScore) {
                                n_better++;
                                // don't export data for matches with low score
                                if (tmpScore >= 20) {
                                    String ref_out = SearchWorker.get_psm_outline(tmpScoreResult);
                                    out_data.add(ref_out);

                                }
                            }

                            if(CParameter.fast_model && n_better >= 1) {
                                break;
                            }
                        }else{
                            if(spectra_matched){
                                // Since the peptide index is already sorted based on peptide mass from small to large,
                                // if there is already spectra matched and the current match has delta mass (del) > CParameter.tol,
                                // we no longer need to compare for the remaining peptides.
                                break;
                            }
                        }
                    }
                }
            }
        }

        // save the best matching from the reference database searching for the target spectrum
        if(bestRefScoreResult.score > 0){
            String ref_out = SearchWorker.get_psm_outline(bestRefScoreResult);
            if(!out_data.contains(ref_out)){
                out_data.add(ref_out);
            }
            if(CParameter.generate_psm_annotation_data) {
                // add PSM annotation
                String bestRef_psmMatch = PeakAnnotation.getPeakAnnotation(bestRefScoreResult.jPeptide.peptide, spectrum, true);
                if (bestRef_psmMatch != null) {
                    psm_annotation_data.add(bestRef_psmMatch);
                }
            }
        }

        this.refPepMatchResult.n_matched_peptides_with_better_score = n_better;
        this.refPepMatchResult.n_matched_peptides = n_matched_peptides;
        this.refPepMatchResult.out_data = out_data;
        this.refPepMatchResult.psm_annotation_data = psm_annotation_data;
    }
}
