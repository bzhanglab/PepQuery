package main.java.OpenModificationSearch;

import main.java.PSMMatch.JPeptide;
import main.java.PSMMatch.ScoreResult;
import main.java.pg.SpectraInput;

import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ConcurrentHashMap;

public class SpectrumMatchWorker implements Runnable{

    public static HashMap<Long, ArrayList<JPeptide>> peptideIndexMap = new HashMap<>();
    public static HashMap<String, Double> spectra2score = new HashMap<>();
    private String spectrum_title;
    private String targeted_peptide;
    private BufferedWriter bWriter;

    public SpectrumMatchWorker(String spectrum_title, String pep, BufferedWriter bWriter){
        this.spectrum_title = spectrum_title;
        this.targeted_peptide = pep;
        this.bWriter = bWriter;

    }

    @Override
    public void run() {
        ConcurrentHashMap<Long, ArrayList<ScoreResult>> res = ModificationDB.getInstance().getCandidateRefPeptidesSingleThread(spectrum_title, peptideIndexMap, targeted_peptide);
        ModificationDB.getInstance().save(bWriter, res, spectra2score, SpectraInput.spectraMap.get(spectrum_title));
    }
}
