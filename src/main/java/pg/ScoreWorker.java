package main.java.pg;

import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import main.java.PSMMatch.JPeptide;
import main.java.PSMMatch.ScoreResult;


import java.util.Iterator;

import static main.java.pg.PeptideSearchMT.scorePeptide2Spectrum;

/**
 * scoring PSM (multiple treads)
 */
public final class ScoreWorker implements  Runnable {

    private PeptideInput peptideInput = null;

    public ScoreWorker(PeptideInput pepInput){
        this.peptideInput = pepInput;
    }



    @Override
    public void run() {
        // save spectra
        for (JPeptide jPeptide : peptideInput.getPtmIsoforms()) {

            Iterator<Spectrum> im = jPeptide.mSnSpectrums.iterator();
            while(im.hasNext()){
                Spectrum spectrum = im.next();
                ScoreResult scoreResult = scorePeptide2Spectrum(jPeptide.peptide,spectrum);

                if(scoreResult.score > CParameter.minScore){
                    jPeptide.scores.add(scoreResult);
                }else{
                    im.remove();
                }

            }

        }

    }
}
