package main.java.OpenModificationSearch;

import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.identification.protein_sequences.digestion.IteratorFactory;
import com.compomics.util.experiment.identification.protein_sequences.digestion.SequenceIterator;
import com.compomics.util.experiment.identification.protein_sequences.digestion.ExtendedPeptide;
import com.compomics.util.parameters.identification.search.DigestionParameters;

import main.java.pg.CParameter;
import main.java.pg.PeptideInput;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.concurrent.ConcurrentHashMap;

public final  class DigestProteinWorker implements Runnable{

    //private ConcurrentHashMap<String,HashMap<String,Double>> pro2pep = null;
    private ConcurrentHashMap<String,ArrayList<Peptide>> pro2pep;
    private ConcurrentHashMap<String,HashSet<String>> pro2pepSeq;
    private String proteinID;
    private String proteinSequence;
    private IteratorFactory iteratorModifications = new IteratorFactory(new ArrayList<>());
    private DigestionParameters digestionParameters;
    private boolean digest_only = false;


    public DigestProteinWorker(String ID, String proSeq, DigestionParameters digestionParameters, ConcurrentHashMap<String,ArrayList<Peptide>> pro2pepMap){
        this.proteinID = ID;
        this.proteinSequence = proSeq;
        this.pro2pep = pro2pepMap;
        this.digestionParameters = digestionParameters;
    }

    public DigestProteinWorker(String ID, String proSeq, DigestionParameters digestionParameters, ConcurrentHashMap<String, HashSet<String>> pro2pepMap,
                               boolean digest_only){
        this.proteinID = ID;
        this.proteinSequence = proSeq;
        this.pro2pepSeq = pro2pepMap;
        this.digestionParameters = digestionParameters;
        this.digest_only = digest_only;
    }

    /**
    public DigestProteinWorker(String ID, String proSeq, ArrayList<String> fixedModifications, DigestionParameters digestionParameters, ConcurrentHashMap<String,ArrayList<Peptide>> pro2pepMap){
        this.proteinID = ID;
        this.proteinSequence = proSeq;
        this.pro2pep = pro2pepMap;
        this.iteratorModifications = new IteratorFactory(fixedModifications);
        this.digestionParameters = digestionParameters;
    }**/

    /**
    public DigestProteinWorker(String ID, String proSeq, DigestionParameters digestionParameters, ConcurrentHashMap<String, HashMap<String, Double>> pro2pepMap){
        this.proteinID = ID;
        this.proteinSequence = proSeq;
        this.pro2pep = pro2pepMap;
        ArrayList<String> fixT = new ArrayList<>();
        this.iteratorModifications = new IteratorFactory(fixT);
        this.digestionParameters = digestionParameters;
    }**/

    @Override
    public void run() {

        proteinSequence = proteinSequence.toUpperCase();

        proteinSequence = proteinSequence.replaceAll("\\*", "");
        // needed
        proteinSequence = proteinSequence.replaceAll("I", "L");


        //HashMap<String,Double> pep2mass = new HashMap<>();
        ArrayList<Peptide> peptides = new ArrayList<>(20);
        HashSet<String> peptideSequences = new HashSet<>(20);
        //SequenceIterator sequenceIterator = iteratorModifications.getSequenceIterator(this.proteinSequence, digestionPreferences, 600.0, 5000.0);
        SequenceIterator sequenceIterator = null;
        try {
            sequenceIterator = iteratorModifications.getSequenceIterator(this.proteinSequence,
                    digestionParameters,
                    CParameter.minPeptideMass,
                    CParameter.maxPeptideMass);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        ExtendedPeptide extendedPeptide;
        while (true) {
            try {
                if ((extendedPeptide = sequenceIterator.getNextPeptide()) == null) {
                    break;
                }

                Peptide peptide = extendedPeptide.peptide;
                // don't use the following way, it's slow.
                //Peptide peptide = new Peptide(extendedPeptide.peptide.getSequence());
                //PeptideInput.addFixedModification(peptide);

                if (peptide.getSequence().length() < CParameter.minPeptideLength || peptide.getSequence().length() > CParameter.maxPeptideLength) {
                    continue;
                }

                if(peptide.getSequence().contains("X")){
                    continue;
                }

                //pep2mass.put(peptide.getSequence(), peptide.getMass());
                if(this.digest_only) {
                    peptideSequences.add(peptide.getSequence());
                }else{
                    ArrayList<Peptide> peptideIsoforms = PeptideInput.calcPeptideIsoforms(peptide.getSequence());
                    peptides.addAll(peptideIsoforms);
                }


                //if (peptide.getSequence().length() >= CParameter.minPeptideLength && peptide.getSequence().length() <= CParameter.maxPeptideLength) {
                    //pep2mass.put(peptide.getSequence(), JPeptide.getMass(peptide));
                    // only sequence and the mass with modification are used.

                    //System.out.println("test:"+peptide.getSequence()+"\t"+JPeptide.getMass(peptide));
                //}

            } catch (InterruptedException e) {
                e.printStackTrace();
            }



        }
        if(this.digest_only){
            pro2pepSeq.put(this.proteinID,peptideSequences);
        }else {
            pro2pep.put(this.proteinID, peptides);
        }

    }
}
