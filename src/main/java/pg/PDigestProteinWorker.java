package main.java.pg;


import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.identification.protein_sequences.digestion.ExtendedPeptide;
import com.compomics.util.experiment.identification.protein_sequences.digestion.IteratorFactory;

import com.compomics.util.experiment.identification.protein_sequences.digestion.SequenceIterator;
import com.compomics.util.parameters.identification.search.DigestionParameters;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.concurrent.ConcurrentHashMap;


/**
 * Digest a protein database. Don't consider modification and only consider enzyme and the max number of cleavage sites
 *
 */
public final  class PDigestProteinWorker implements Runnable{

    private ConcurrentHashMap<String,HashSet<String>> pro2pep = null;
    private String proteinID;
    private String proteinSequence;
    private IteratorFactory iteratorModifications;
    private DigestionParameters digestionParameters;
    private ArrayList<String> fixedModifications = new ArrayList<>();

    public PDigestProteinWorker(String ID, String proSeq, DigestionParameters digestionPreferences, ConcurrentHashMap<String, HashSet<String>> pro2pepMap){
        this.proteinID = ID;
        this.proteinSequence = proSeq;
        this.pro2pep = pro2pepMap;
        this.iteratorModifications = new IteratorFactory(fixedModifications);
        this.digestionParameters = digestionPreferences;
    }

    @Override
    public void run() {

        proteinSequence = proteinSequence.toUpperCase();

        proteinSequence = proteinSequence.replaceAll("\\*", "");
        proteinSequence = proteinSequence.replaceAll("I", "L");


        HashSet<String> pepSet = new HashSet<>();

        // SequenceIterator sequenceIterator = iteratorModifications.getSequenceIterator(this.proteinSequence, digestionPreferences, 600.0, 5000.0);
        SequenceIterator sequenceIterator = null;
        try {
            sequenceIterator = iteratorModifications.getSequenceIterator(this.proteinSequence, digestionParameters, CParameter.minPeptideMass, CParameter.maxPeptideMass);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        ExtendedPeptide peptideWithPosition;
        while (true) {
            try {
                if (!((peptideWithPosition = sequenceIterator.getNextPeptide()) != null)) {
                    break;
                }

                Peptide peptide = peptideWithPosition.peptide;

                if (peptide.getSequence().length() >= CParameter.minPeptideLength && peptide.getSequence().length() <= CParameter.maxPeptideLength) {
                    pepSet.add(peptide.getSequence());
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            }



        }

        pro2pep.put(this.proteinID,pepSet);

    }
}
