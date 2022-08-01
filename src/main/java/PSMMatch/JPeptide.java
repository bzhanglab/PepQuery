package main.java.PSMMatch;

import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import com.compomics.util.parameters.identification.search.ModificationParameters;
import main.java.pg.JSequenceProvider;

import java.util.ArrayList;
import java.util.HashMap;


public class JPeptide {

    public Peptide peptide;
    private double mass = -1.0;
    public double isotope_error = 0;

    public ArrayList<Spectrum> mSnSpectrums = new ArrayList<>();
    public ArrayList<String> spectraIndexs = new ArrayList<>();

    /**
     * Because sometimes there is no charge information for a spectrum in mgf file, here we need to save the charge
     * information when perform matching.
     */
    // public ArrayList<Integer> charges = new ArrayList<>(1);
    public HashMap<String,Integer> spectrumID2charge = new HashMap<>();

    public ArrayList<ScoreResult> scores = new ArrayList<>();

    // Whether the match is valid. If it's not valid, then it will not be used for downstream processing.
    public ArrayList<Boolean> valid = new ArrayList<>();

    public JPeptide(Peptide pep){
        this.peptide = pep;
    }


    public synchronized void saveMatchedSpectrum(ScoreResult scoreResult,String spectrumID, boolean isValid){
        this.scores.add(scoreResult);
        this.spectraIndexs.add(spectrumID);
        this.valid.add(isValid);
    }

    public synchronized void saveMatchedSpectrum(ScoreResult scoreResult,String spectrumID, int charge, boolean isValid){
        this.scores.add(scoreResult);
        this.spectraIndexs.add(spectrumID);
        spectrumID2charge.put(spectrumID,charge);
        this.valid.add(isValid);
    }


    public double getMass(){
        if(mass > 0){
            return mass;
        }else{
            mass = getMass(this.peptide);
            return mass;
        }

    }


    private static final ModificationParameters modificationParameters = new ModificationParameters();
    private static final JSequenceProvider sequenceProvider = new JSequenceProvider();
    private static final SequenceMatchingParameters sequenceMatchingParameters = new SequenceMatchingParameters();

    /**
     * Use this function to get peptide mass not the peptide.getMass
     * @param pep
     * @return
     */
    public static double getMass(Peptide pep){
        // utilities uses -1.0.
        pep.setMass(-1.0);
        double pep_mass = pep.getMass(modificationParameters,sequenceProvider,sequenceMatchingParameters);
        return pep_mass;
    }


}
