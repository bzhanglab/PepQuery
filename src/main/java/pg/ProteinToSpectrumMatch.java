package main.java.pg;

import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.identification.protein_sequences.digestion.ExtendedPeptide;
import com.compomics.util.experiment.identification.protein_sequences.digestion.IteratorFactory;
import com.compomics.util.experiment.identification.protein_sequences.digestion.SequenceIterator;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import com.compomics.util.parameters.identification.search.DigestionParameters;
import main.java.PSMMatch.JPeptide;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ConcurrentHashMap;

public class ProteinToSpectrumMatch implements  Runnable{

    private IteratorFactory iteratorFactory;
    private String proSeq;
    private DigestionParameters digestionParameters;

    // if a peptide is searched before, then it will be omitted.
    public ConcurrentHashMap<String,Integer> searchedPeptides = new ConcurrentHashMap<>();

    // prepare spectrum index for fast peptide query
    HashMap<Integer,ArrayList<String>> spectraIndex;

    /**
     * key: spectrum title, value: an ArrayList of Peptide objects
     */
    public HashMap<String,ArrayList<Peptide>> ms2peptide;


    public ProteinToSpectrumMatch (IteratorFactory iteratorFactory,
                                   String proSeq,
                                   DigestionParameters digestionParameters,
                                   ConcurrentHashMap<String,Integer> searchedPeptides,
                                   HashMap<Integer,ArrayList<String>> spectraIndex,
                                   HashMap<String,ArrayList<Peptide>> ms2peptide){
        this.iteratorFactory = iteratorFactory;
        this.proSeq = proSeq.replaceAll("I","L");
        this.digestionParameters = digestionParameters;
        this.searchedPeptides = searchedPeptides;
        this.spectraIndex = spectraIndex;
        this.ms2peptide = ms2peptide;
    }


    @Override
    public void run() {

        SequenceIterator sequenceIterator = null;
        try {
            sequenceIterator = iteratorFactory.getSequenceIterator(proSeq, digestionParameters, CParameter.minPeptideMass, CParameter.maxPeptideMass);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        int num = 0;
        ExtendedPeptide peptideWithPosition = new ExtendedPeptide();

        Spectrum spectrum;

        double peptideMass;
        int va;
        double mass;
        double del;

        while (true) {
            try {
                if ((peptideWithPosition = sequenceIterator.getNextPeptide()) == null) break;
            } catch (InterruptedException e) {
                e.printStackTrace();
            }

            Peptide peptide = peptideWithPosition.peptide;


            if (peptide.getSequence().length() < CParameter.minPeptideLength || peptide.getSequence().length() > CParameter.maxPeptideLength) {
                continue;
            }

            if(searchedPeptides.containsKey(peptide.getSequence())){
                continue;
            }else{
                searchedPeptides.put(peptide.getSequence(),1);
            }

            if(peptide.getSequence().contains("X")){
                continue;
            }

            // Consider modification
            ArrayList<Peptide> peptideIsoforms = PeptideInput.calcPeptideIsoforms(peptide.getSequence());
            //peptideIsoforms.add(peptide);
            for(Peptide p: peptideIsoforms){
                JPeptide jPeptide = new JPeptide(p);
                peptideMass = jPeptide.getMass();
                // consider isotope error
                ArrayList<Double> isotope_masses = CParameter.get_peptide_mass_with_isotope_error(peptideMass);
                for(double i_mass : isotope_masses) {
                    va = (int) Math.round(10.0 * i_mass);
                    // +/- 0.1 Da
                    for (int i = (va - 1); i <= (va + 1); i++) {
                        if (spectraIndex.containsKey(i)) {
                            for (String title : spectraIndex.get(i)) {
                                spectrum = SpectraInput.spectraMap.get(title);
                                mass = spectrum.getPrecursor().getMass(spectrum.getPrecursor().possibleCharges[0]);
                                // use peptide mass not i_mass.
                                ArrayList<Double> tol_res = CParameter.get_mass_error(mass,peptideMass);
                                if(tol_res.size()>=1) {
                                    //del = Math.abs(i_mass - mass);
                                    if (CParameter.tolu.equalsIgnoreCase("ppm")) {
                                        del = tol_res.get(0);
                                    }else{
                                        del = tol_res.get(1);
                                    }
                                    if (del <= CParameter.tol) {
                                        //System.out.println("spectra match:"+title+"\t"+mass+"\t"+peptideMass+"\t"+p.getSequence()+"\t"+ ModificationDB.getModificationString(p));
                                        if (ms2peptide.containsKey(title)) {
                                            ms2peptide.get(title).add(p);
                                        } else {
                                            ArrayList<Peptide> peps = new ArrayList<>();
                                            peps.add(p);
                                            ms2peptide.put(title, peps);
                                        }
                                    }
                                }

                            }
                        }
                    }
                }
            }

        }


        if( (num%1000)==0){
            //Cloger.getInstance().logger.info("Finished proteins:"+num);

        }
    }
}
