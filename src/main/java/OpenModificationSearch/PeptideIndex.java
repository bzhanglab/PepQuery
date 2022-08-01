package main.java.OpenModificationSearch;

import com.compomics.util.experiment.biology.proteins.Peptide;

import java.io.Serializable;
import java.util.ArrayList;


public class PeptideIndex implements Serializable {


    private static final long serialVersionUID = 123L;

    public ArrayList<String> peptideSequences = new ArrayList<>();
    public ArrayList<Double> peptideMasses = new ArrayList<>();
    public ArrayList<Peptide> peptides = new ArrayList<>();



    public PeptideIndex(String seq, double mass){
        this.peptideSequences.add(seq);
        this.peptideMasses.add(mass);
    }

    public PeptideIndex(Peptide pep, double mass){
        this.peptides.add(pep);
        this.peptideMasses.add(mass);
    }

    public PeptideIndex(Peptide pep){
        this.peptides.add(pep);
    }

}
