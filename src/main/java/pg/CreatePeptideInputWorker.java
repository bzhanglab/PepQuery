package main.java.pg;


import com.compomics.util.experiment.biology.proteins.Peptide;


import java.util.ArrayList;
import java.util.concurrent.ConcurrentHashMap;

public class CreatePeptideInputWorker implements  Runnable {

    private String peptideSequence;

    public ConcurrentHashMap<String,PeptideInput> peptideInputConcurrentHashMap = new ConcurrentHashMap<>();

    public CreatePeptideInputWorker(String pep, ConcurrentHashMap<String,PeptideInput> pMap){
        this.peptideSequence = pep;
        this.peptideInputConcurrentHashMap = pMap;
    }

    @Override
    public void run() {
        PeptideInput peptideInput = new PeptideInput(this.peptideSequence);
        ArrayList<Peptide> peptides = PeptideInput.calcPeptideIsoforms(peptideSequence);
        peptideInput.setPtmIsoforms(peptides);
        peptideInputConcurrentHashMap.put(peptideSequence,peptideInput);

    }
}
