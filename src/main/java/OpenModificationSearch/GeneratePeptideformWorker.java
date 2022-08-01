package main.java.OpenModificationSearch;

import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import main.java.PSMMatch.JPeptide;
import main.java.pg.PeptideInput;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ConcurrentHashMap;

public class GeneratePeptideformWorker implements Runnable{

    private static double min_mass = 0.0;
    private static double max_mass = 100000.0;

    private String peptide_sequence;
    private ConcurrentHashMap<String, ArrayList<JPeptide>> searchedPeptides;
    public GeneratePeptideformWorker(String seq, ConcurrentHashMap<String, ArrayList<JPeptide>> searchedPeptides){
        this.peptide_sequence = seq;
        this.searchedPeptides = searchedPeptides;
    }

    @Override
    public void run() {
        ArrayList<Peptide> peptideIsoforms = PeptideInput.calcPeptideIsoforms(this.peptide_sequence);
        for(Peptide peptide:peptideIsoforms){
            JPeptide jPeptide = new JPeptide(peptide);
            double mass = jPeptide.getMass();
            if(mass > min_mass && mass < max_mass) {
                this.searchedPeptides.get(this.peptide_sequence).add(jPeptide);
            }
        }

    }

    public static void reset_mass_range(){
        min_mass = 0.0;
        max_mass = 100000.0;
    }

    public static void set_mass_range(HashMap<String, Spectrum> spectraMap){
        double min_spectrum_mass = Double.MAX_VALUE;
        double max_spectrum_mass = 0;
        for(Spectrum spectrum:spectraMap.values()){
            double mass = spectrum.getPrecursor().getMass(spectrum.getPrecursor().possibleCharges[0]);
            min_spectrum_mass = Math.min(min_spectrum_mass,mass);
            max_spectrum_mass = Math.max(max_spectrum_mass,mass);
        }
        min_mass = min_spectrum_mass - Math.abs(ModificationDB.leftMass) - 5.0;
        max_mass = max_spectrum_mass + Math.abs(ModificationDB.rightMass) + 5.0;
    }
}
