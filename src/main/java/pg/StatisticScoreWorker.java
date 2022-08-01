package main.java.pg;


import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import main.java.OpenModificationSearch.ModificationDB;
import main.java.plot.PeakAnnotation;
import main.java.PSMMatch.JPeptide;
import main.java.PSMMatch.ScoreResult;
import main.java.util.Cloger;
import java.util.ArrayList;

import static main.java.pg.PeptideSearchMT.scorePeptide2Spectrum;

/**
 * Scoring: single spectrum and single peptide.
 * Speed up the whole processing when taking as input of single peptide for identification.
 */
public class StatisticScoreWorker implements  Runnable{

    private PeptideInput peptideInput = null;
    private int i=0;

    private String fixModOptionValue = null;

    private JPeptide jPeptide = null;

    public StatisticScoreWorker (PeptideInput pi, JPeptide jPeptide, int index){
        this.peptideInput = pi;
        this.i = index;
        this.jPeptide = jPeptide;

    }


    @Override
    public void run() {
        Cloger.getInstance().logger.debug(Thread.currentThread().getName()+" -> Start "+this.i+" ...");

        //MSnSpectrum spectrum = jPeptide.mSnSpectrums.get(i);
        Spectrum spectrum = SpectraInput.spectraMap.get(jPeptide.spectraIndexs.get(i));
        ScoreResult scoreResult = jPeptide.scores.get(i);
        if(peptideInput.randomPeptides.size()>=1) {
            int nRand = 0;

            ScoreResult tmpScoreResult = new ScoreResult();
            for (String randomPeptideSeq : peptideInput.randomPeptides) {
                Peptide randomPeptide = PeptideInput.getModificationPeptide(jPeptide.peptide, randomPeptideSeq, fixModOptionValue);
                tmpScoreResult = scorePeptide2Spectrum(randomPeptide, spectrum);
                double randomScore = tmpScoreResult.score;
                if (randomScore >= scoreResult.score) {
                    nRand++;
                }

            }
            scoreResult.p_value = 1.0 * (nRand + 1) / (peptideInput.randomPeptides.size() + 1);
            scoreResult.n_total_random_peptides = peptideInput.randomPeptides.size();
            scoreResult.n_random_better_match = nRand;
        }else{
            Cloger.getInstance().logger.error("There is no random peptide!");
            System.exit(1);
        }
    }
}
