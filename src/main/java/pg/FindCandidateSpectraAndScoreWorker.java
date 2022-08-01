package main.java.pg;

import com.compomics.util.experiment.mass_spectrometry.spectra.Precursor;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import main.java.OpenModificationSearch.ModificationDB;
import main.java.PSMMatch.JPeptide;
import main.java.PSMMatch.ScoreResult;
import main.java.util.Cloger;

import java.util.ArrayList;
import java.util.concurrent.ConcurrentHashMap;
import static main.java.pg.PeptideSearchMT.scorePeptide2Spectrum;

/**
 * Search a spectrum against peptides (peptide index). Get the matched peptides and score PSMs. Set a min score cutoff
 * for matching.
 */
public class FindCandidateSpectraAndScoreWorker implements  Runnable {


    public ConcurrentHashMap<Integer,ArrayList<JPeptide>> pIndex = new ConcurrentHashMap<>();
    public Spectrum spectrum = new Spectrum();


    public FindCandidateSpectraAndScoreWorker(ConcurrentHashMap<Integer,ArrayList<JPeptide>> pepIndex, Spectrum msnspectrum){
        this.pIndex = pepIndex;
        this.spectrum = msnspectrum;

    }


    @Override
    public void run() {

        // System.out.println("Spectrum:"+spectrum.getSpectrumTitle());

        if(spectrum.getPrecursor().possibleCharges.length>=1) {

            int charge = spectrum.getPrecursor().possibleCharges[0];
            double exp_mass = spectrum.getPrecursor().getMass(charge);

            ArrayList<Double> isotope_masses = CParameter.get_precursor_mass_with_isotope_error(exp_mass);

            for(double mass:isotope_masses) {
                //Cloger.getInstance().logger.info(exp_mass+"\t"+mass);

                int va = (int) Math.round(mass * 10.0);
                double peptideMass;
                double del;
                for (int i = (va - 1); i <= (va + 1); i++) {
                    if (pIndex.containsKey(i)) {
                        for (JPeptide jPeptide : pIndex.get(i)) {
                            peptideMass = jPeptide.peptide.getMass();
                            del = Math.abs(peptideMass - mass);
                            if (CParameter.tolu.equalsIgnoreCase("ppm")) {
                                del = (1.0e6) * 1.0 * del / peptideMass;
                            }
                            if (del <= CParameter.tol) {

                                if(SpectraInput.pep2targeted_spectra.size()>=1 &&
                                        SpectraInput.pep2targeted_spectra.containsKey(jPeptide.peptide.getSequence()) &&
                                        !SpectraInput.pep2targeted_spectra.get(jPeptide.peptide.getSequence()).containsKey(spectrum.getSpectrumTitle())){
                                    // If pep2targeted_spectra is not empty, it means the input (-pep) also contains spectrum title
                                    // If for the target peptide, the current spectrum (spectrum object) is not present in
                                    // pep2targeted_spectra, this spectrum will be ignored.
                                    continue;
                                }

                                // If SpectraInput.pep2targeted_spectra contains the current peptide and the spectrum, we need
                                // to update the status of the targeted PSM
                                if(SpectraInput.pep2targeted_spectra.containsKey(jPeptide.peptide.getSequence()) &&
                                        SpectraInput.pep2targeted_spectra.get(jPeptide.peptide.getSequence()).containsKey(spectrum.getSpectrumTitle())){
                                    SpectraInput.pep2targeted_spectra.get(jPeptide.peptide.getSequence()).get(spectrum.getSpectrumTitle()).update_status(PsmType.spectrum_matched);
                                }

                                ScoreResult scoreResult = scorePeptide2Spectrum(jPeptide.peptide, spectrum);

                                if (scoreResult.score > CParameter.minScore) {

                                    // make sure there is only one thread can access it
                                    jPeptide.saveMatchedSpectrum(scoreResult, spectrum.getSpectrumTitle(), true);
                                    //jPeptide.scores.add(scoreResult);
                                    //jPeptide.spectraIndexs.add(spectrum.getSpectrumTitle());
                                    //jPeptide.valid.add(true);
                                    //jPeptide.mSnSpectrums.add(spectrum);
                                }else{
                                    if(SpectraInput.pep2targeted_spectra.containsKey(jPeptide.peptide.getSequence()) &&
                                            SpectraInput.pep2targeted_spectra.get(jPeptide.peptide.getSequence()).containsKey(spectrum.getSpectrumTitle())){
                                        Cloger.getInstance().logger.warn("Low score for targeted PSM: "+
                                                jPeptide.peptide.getSequence()+"\t"+
                                                ModificationDB.getInstance().getModificationString(jPeptide.peptide)+"\t"+
                                                spectrum.getSpectrumTitle()+"\t"+
                                                scoreResult.score);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }else{

            // if no charge state for this spectrum
            int charge = 2;

            // We need to change the precursor of this spectrum

            int possibleCharges[] = new int[]{charge};
            Precursor precursor = new Precursor(0.0,spectrum.getPrecursor().mz,possibleCharges);
            spectrum.precursor = precursor;

            double exp_mass = spectrum.getPrecursor().getMass(charge);

            ArrayList<Double> isotope_masses = CParameter.get_precursor_mass_with_isotope_error(exp_mass);

            boolean find = false;
            for(double mass:isotope_masses) {
                //Cloger.getInstance().logger.info(exp_mass+"\t"+mass);
                int va = (int) Math.round(mass * 10.0);
                for (int i = (va - 1); i <= (va + 1); i++) {
                    if (pIndex.containsKey(i)) {
                        for (JPeptide jPeptide : pIndex.get(i)) {
                            double peptideMass = jPeptide.peptide.getMass();
                            double del = Math.abs(peptideMass - mass);
                            if (CParameter.tolu.equalsIgnoreCase("ppm")) {
                                del = (1.0e6) * 1.0 * Math.abs(peptideMass - mass) / peptideMass;
                            }
                            if (del <= CParameter.tol) {

                                if(SpectraInput.pep2targeted_spectra.size()>=1 &&
                                        SpectraInput.pep2targeted_spectra.containsKey(jPeptide.peptide.getSequence()) &&
                                        !SpectraInput.pep2targeted_spectra.get(jPeptide.peptide.getSequence()).containsKey(spectrum.getSpectrumTitle())){
                                    // If pep2targeted_spectra is not empty, it means the input (-pep) also contains spectrum title
                                    // If for the target peptide, the current spectrum (spectrum object) is not present in
                                    // pep2targeted_spectra, this spectrum will be ignored.
                                    continue;
                                }

                                ScoreResult scoreResult = scorePeptide2Spectrum(jPeptide.peptide, spectrum);

                                if (scoreResult.score > CParameter.minScore) {

                                    // jPeptide.saveMatchedSpectrum(scoreResult,spectrum.getSpectrumTitle(),true);
                                    jPeptide.saveMatchedSpectrum(scoreResult, spectrum.getSpectrumTitle(), charge, true);
                                    //jPeptide.scores.add(scoreResult);
                                    //jPeptide.spectraIndexs.add(spectrum.getSpectrumTitle());
                                    //jPeptide.valid.add(true);
                                    //jPeptide.mSnSpectrums.add(spectrum);
                                    find = true;
                                }else{
                                    if(SpectraInput.pep2targeted_spectra.containsKey(jPeptide.peptide.getSequence()) &&
                                            SpectraInput.pep2targeted_spectra.get(jPeptide.peptide.getSequence()).containsKey(spectrum.getSpectrumTitle())){
                                        Cloger.getInstance().logger.warn("Low score for targeted PSM: "+
                                                jPeptide.peptide.getSequence()+"\t"+
                                                ModificationDB.getInstance().getModificationString(jPeptide.peptide)+"\t"+
                                                spectrum.getSpectrumTitle()+"\t"+
                                                scoreResult.score);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // If no result for charge state 2, will consider 3
            if(!find) {
                charge = 3;

                // We need to change the precursor of this spectrum
                possibleCharges = new int[]{charge};
                Precursor precursor3 = new Precursor(0.0,spectrum.getPrecursor().mz,possibleCharges);
                spectrum.precursor = precursor3;

                exp_mass = spectrum.getPrecursor().getMass(charge);
                isotope_masses = CParameter.get_precursor_mass_with_isotope_error(exp_mass);
                for(double mass:isotope_masses) {
                    //Cloger.getInstance().logger.info(exp_mass+"\t"+mass);
                    int va = (int) Math.round(mass * 10.0);

                    for (int i = (va - 1); i <= (va + 1); i++) {
                        if (pIndex.containsKey(i)) {
                            for (JPeptide jPeptide : pIndex.get(i)) {
                                double peptideMass = jPeptide.peptide.getMass();
                                double del = Math.abs(peptideMass - mass);
                                if (CParameter.tolu.equalsIgnoreCase("ppm")) {
                                    del = (1.0e6) * 1.0 * Math.abs(peptideMass - mass) / peptideMass;
                                }
                                if (del <= CParameter.tol) {

                                    if(SpectraInput.pep2targeted_spectra.size()>=1 &&
                                            SpectraInput.pep2targeted_spectra.containsKey(jPeptide.peptide.getSequence()) &&
                                            !SpectraInput.pep2targeted_spectra.get(jPeptide.peptide.getSequence()).containsKey(spectrum.getSpectrumTitle())){
                                        // If pep2targeted_spectra is not empty, it means the input (-pep) also contains spectrum title
                                        // If for the target peptide, the current spectrum (spectrum object) is not present in
                                        // pep2targeted_spectra, this spectrum will be ignored.
                                        continue;
                                    }

                                    ScoreResult scoreResult = scorePeptide2Spectrum(jPeptide.peptide, spectrum);

                                    if (scoreResult.score > CParameter.minScore) {

                                        // jPeptide.saveMatchedSpectrum(scoreResult,spectrum.getSpectrumTitle(),true);
                                        jPeptide.saveMatchedSpectrum(scoreResult, spectrum.getSpectrumTitle(), charge, true);
                                        //jPeptide.scores.add(scoreResult);
                                        //jPeptide.spectraIndexs.add(spectrum.getSpectrumTitle());
                                        //jPeptide.valid.add(true);
                                        //jPeptide.mSnSpectrums.add(spectrum);
                                    }else{
                                        if(SpectraInput.pep2targeted_spectra.containsKey(jPeptide.peptide.getSequence()) &&
                                                SpectraInput.pep2targeted_spectra.get(jPeptide.peptide.getSequence()).containsKey(spectrum.getSpectrumTitle())){
                                            Cloger.getInstance().logger.warn("Low score for targeted PSM: "+
                                                    jPeptide.peptide.getSequence()+"\t"+
                                                    ModificationDB.getInstance().getModificationString(jPeptide.peptide)+"\t"+
                                                    spectrum.getSpectrumTitle()+"\t"+
                                                    scoreResult.score);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

            }

        }

    }
}
