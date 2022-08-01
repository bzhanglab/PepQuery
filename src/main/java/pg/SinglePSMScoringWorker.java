package main.java.pg;

import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import main.java.PSMMatch.JPeptide;
import main.java.PSMMatch.ScoreResult;

import java.lang.reflect.Array;
import java.util.ArrayList;

import static main.java.pg.PeptideSearchMT.scorePeptide2Spectrum;

public class SinglePSMScoringWorker implements Runnable {


    private ArrayList<JPeptide> peptides = new ArrayList<>();
    private Spectrum spectrum;

    public SinglePSMScoringWorker(ArrayList<JPeptide> peptides, Spectrum spectrum){
        this.peptides = peptides;
        this.spectrum = spectrum;
    }

    @Override
    public void run() {
        for(JPeptide peptide : this.peptides){
            ScoreResult scoreResult = scorePeptide2Spectrum(peptide.peptide, spectrum);
        }
        // cannot use this because the same peptide object could be used by different spectra
        // scores is not thread safe
        //this.peptide.scores.add(scoreResult);
    }
}
