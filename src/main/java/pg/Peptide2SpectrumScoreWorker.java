package main.java.pg;


import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import main.java.PSMMatch.JPeptide;
import main.java.PSMMatch.ScoreResult;


import java.util.concurrent.ConcurrentHashMap;

import static main.java.pg.PeptideSearchMT.scorePeptide2Spectrum;

/**
 * When input MS/MS data is in SQL format. Search target peptide against the MS/MS data and score PSMs.
 */
public class Peptide2SpectrumScoreWorker implements  Runnable{

    public Spectrum spectrum = new Spectrum();
    public JPeptide jPeptide;

    /**
     * Save spectra when using multiple threads
     */
    ConcurrentHashMap<String,Spectrum> matchedSpectra = new ConcurrentHashMap<>();


    public Peptide2SpectrumScoreWorker(JPeptide peptide, Spectrum msnspectrum,ConcurrentHashMap<String,Spectrum> mSpectra){
        this.jPeptide = peptide;
        this.spectrum = msnspectrum;
        this.matchedSpectra = mSpectra;

    }

    @Override
    public void run() {
        ScoreResult scoreResult = scorePeptide2Spectrum(jPeptide.peptide,spectrum);

        if(scoreResult.score > CParameter.minScore){
            // save scoring result
            jPeptide.scores.add(scoreResult);
            jPeptide.spectraIndexs.add(spectrum.getSpectrumTitle());
            jPeptide.valid.add(true);
            //jPeptide.mSnSpectrums.add(spectrum);

            if(!matchedSpectra.containsKey(spectrum.getSpectrumTitle())){
                matchedSpectra.put(spectrum.getSpectrumTitle(),spectrum);
            }
        }
    }
}
