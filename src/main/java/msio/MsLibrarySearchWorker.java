package main.java.msio;


import com.compomics.util.experiment.io.mass_spectrometry.mgf.MgfFileIterator;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import com.compomics.util.gui.waiting.waitinghandlers.WaitingHandlerCLIImpl;
import com.compomics.util.waiting.WaitingHandler;
import main.java.OpenModificationSearch.ModificationDB;
import main.java.PSMMatch.ScoreResult;
import main.java.pg.*;
import main.java.util.Cloger;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;


import static main.java.index.IndexWorker.decompressGzipFile;
import static main.java.index.IndexWorker.decompressXzipFile;
import static main.java.pg.PeptideSearchMT.scorePeptide2Spectrum;
import static main.java.pg.SpectraInput.getRangeOfMass;

/**
 * For each target peptide, find the candidate spectra and perform scoring.
 * Consider different forms (different potential modification combinations) of the peptide sequence
 */
public class MsLibrarySearchWorker implements Runnable {

    private PeptideInput peptideInput;
    private String msms_library_path;

    public MsLibrarySearchWorker(PeptideInput pepInput, String msms_path){
        this.peptideInput = pepInput;
        this.msms_library_path = msms_path;
        if(this.msms_library_path.startsWith("s3:")){
            this.msms_library_path = this.msms_library_path.replaceFirst("/$","");
        }
    }

    @Override
    public void run() {

        // key: peptide mass, value: peptide form in peptideInput.getPtmIsoforms()
        // This is used to speed up spectra finding. For the peptide forms with same mass, only one time finding is
        // needed.
        HashMap<Double, ArrayList<Integer>> pepMass2index = new HashMap<>();
        for (int j = 0; j < peptideInput.getPtmIsoforms().size(); j++) {
            // for each peptide form, get spectra
            ArrayList<Spectrum> rsSpectra = new ArrayList<>();
            double pepIsoformMass = peptideInput.getPtmIsoforms().get(j).getMass();
            if(pepMass2index.containsKey(pepIsoformMass)){
                pepMass2index.get(pepIsoformMass).add(j);
            }else{
                ArrayList<Integer> indexList = new ArrayList<>();
                indexList.add(j);
                pepMass2index.put(pepIsoformMass,indexList);
            }
        }

        //for (int j = 0; j < peptideInput.getPtmIsoforms().size(); j++) {
        for(double pepIsoformMass: pepMass2index.keySet()){
            // for each unique peptide mass, get spectra
            ArrayList<Spectrum> rsSpectra = new ArrayList<>();
            // double pepIsoformMass = peptideInput.getPtmIsoforms().get(j).getMass();
            ArrayList<Double> isotope_masses = CParameter.get_peptide_mass_with_isotope_error(pepIsoformMass);
            for(double i_mass : isotope_masses) {


                double massRange[] = getRangeOfMass(i_mass, CParameter.tol, CParameter.tolu.equalsIgnoreCase("ppm"));
                int left_range = (int) Math.round(massRange[0] * 10);
                int right_range = (int) Math.round(massRange[1] * 10);

                //System.out.println(i_mass+"\t"+massRange[0]+"\t"+massRange[1]);

                for (int i = left_range; i <= right_range; i++) {
                    String i_mgf = msms_library_path + File.separator + i + ".mgf";

                    if (this.msms_library_path.startsWith("s3:")) {
                        String objPath = msms_library_path + "/" + i + ".mgf";
                        File tempFile = null;
                        try {
                            tempFile = File.createTempFile(Thread.currentThread().getName(), "_" + i + ".mgf");
                        } catch (IOException e) {
                            e.printStackTrace();
                            System.exit(1);
                        }

                        i_mgf = tempFile.getAbsolutePath();

                        if (S3Interface.getInstance().download(objPath, i_mgf)) {
                            rsSpectra.addAll(loadSpectra(i_mgf, massRange, CParameter.minPeaks));
                            tempFile.delete();
                        } else {
                            Cloger.getInstance().logger.warn(objPath + " doesn't exist");
                            tempFile.delete();
                        }

                    } else {
                        rsSpectra.addAll(loadSpectra(i_mgf, massRange, CParameter.minPeaks));
                    }
                }
            }

            for(Spectrum spectrum: rsSpectra){

                if(SpectraInput.pep2targeted_spectra.size()>=1 &&
                        SpectraInput.pep2targeted_spectra.containsKey(peptideInput.peptideSequence) &&
                        !SpectraInput.pep2targeted_spectra.get(peptideInput.peptideSequence).containsKey(spectrum.getSpectrumTitle())){
                    // If pep2targeted_spectra is not empty, it means the input (-pep) also contains spectrum title
                    // If for the target peptide, the current spectrum (spectrum object) is not present in
                    // pep2targeted_spectra, this spectrum will be ignored.
                    continue;
                }

                for(int j: pepMass2index.get(pepIsoformMass)) {

                    ArrayList<Double> tol_res = CParameter.get_mass_error(spectrum.getPrecursor().getMass(spectrum.getPrecursor().possibleCharges[0]), peptideInput.getPtmIsoforms().get(j).peptide.getMass());
                    if (tol_res.size() >= 1) {
                        ScoreResult scoreResult = scorePeptide2Spectrum(peptideInput.getPtmIsoforms().get(j).peptide, spectrum);

                        // If SpectraInput.pep2targeted_spectra contains the current peptide and the spectrum, we need
                        // to update the status of the targeted PSM
                        if(SpectraInput.pep2targeted_spectra.containsKey(peptideInput.getPtmIsoforms().get(j).peptide.getSequence()) &&
                                SpectraInput.pep2targeted_spectra.get(peptideInput.getPtmIsoforms().get(j).peptide.getSequence()).containsKey(spectrum.getSpectrumTitle())){
                            SpectraInput.pep2targeted_spectra.get(peptideInput.getPtmIsoforms().get(j).peptide.getSequence()).get(spectrum.getSpectrumTitle()).update_status(PsmType.spectrum_matched);
                        }

                        if (scoreResult.score > CParameter.minScore) {
                            // save scoring result
                            peptideInput.getPtmIsoforms().get(j).scores.add(scoreResult);
                            peptideInput.getPtmIsoforms().get(j).spectraIndexs.add(spectrum.getSpectrumTitle());
                            peptideInput.getPtmIsoforms().get(j).valid.add(true);
                            //jPeptide.mSnSpectrums.add(spectrum);

                            peptideInput.addSpectrum(spectrum);
                        } else {
                            if (SpectraInput.pep2targeted_spectra.containsKey(peptideInput.getPtmIsoforms().get(j).peptide.getSequence()) &&
                                    SpectraInput.pep2targeted_spectra.get(peptideInput.getPtmIsoforms().get(j).peptide.getSequence()).containsKey(spectrum.getSpectrumTitle())) {
                                Cloger.getInstance().logger.warn("Low score for targeted PSM: " +
                                        peptideInput.getPtmIsoforms().get(j).peptide.getSequence() + "\t" +
                                        ModificationDB.getInstance().getModificationString(peptideInput.getPtmIsoforms().get(j).peptide) + "\t" +
                                        spectrum.getSpectrumTitle() + "\t" +
                                        scoreResult.score);
                            }
                        }
                    }
                }

            }

        }
    }

    public static ArrayList<Spectrum> loadSpectra(String mgf, double[] massRange, int minPeaks){

        //System.out.println(mgf);
        ArrayList<Spectrum> spectrums = new ArrayList<>();
        String ms_file = mgf;
        File F = new File(ms_file);
        boolean isZip = false;
        if(!F.isFile()){
            // GZ file
            ms_file = mgf + ".gz";
            F = new File(ms_file);
            if(F.isFile()){
                isZip = true;
                File rawF = new File(ms_file);
                File tempFile = null;
                try {
                    tempFile = File.createTempFile(Thread.currentThread().getName(), rawF.getName().replaceAll(".gz$", ""));
                } catch (IOException e) {
                    e.printStackTrace();
                }
                try {
                    decompressGzipFile(ms_file, tempFile.getAbsolutePath());
                } catch (IOException e) {
                    e.printStackTrace();
                }
                ms_file = tempFile.getAbsolutePath();
            }else{
                ms_file = mgf + ".xz";
                F = new File(ms_file);
                if(F.isFile()) {
                    isZip = true;
                    File rawF = new File(ms_file);
                    File tempFile = null;
                    try {
                        tempFile = File.createTempFile(Thread.currentThread().getName(), rawF.getName().replaceAll(".xz$", ""));
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    try {
                        decompressXzipFile(ms_file, tempFile.getAbsolutePath());
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    ms_file = tempFile.getAbsolutePath();
                }else{
                    // no matched mgf file
                    Cloger.getInstance().logger.warn("File doesn't exist:"+mgf);
                    return(spectrums);
                }
            }
        }

        F = new File(ms_file);
        WaitingHandler waitingHandler = new WaitingHandlerCLIImpl();
        waitingHandler.setDisplayProgress(false);

        MgfFileIterator mgfFileIterator = new MgfFileIterator(F, waitingHandler);
        String title;
        Spectrum spectrum = new Spectrum();
        double mass;
        int charge;
        long spectrum_index = 0;
        String spectrumTitle;
        while ((title = mgfFileIterator.next()) != null) {

            spectrum = mgfFileIterator.getSpectrum();

            // This is for the mgf file from indexed MS/MS spectra library, so each spectrum in the mgf file should have
            // a valid and unique spectrum title.

            // If SpectraInput.spectra2targeted_peptide is not empty, it means the input for -pep contains spectra for
            // each targeted peptide. If this is the case, we only need to score peptides against the targeted spectra
            if(SpectraInput.spectra2targeted_peptide.size()>=1 &&
                    !SpectraInput.spectra2targeted_peptide.containsKey(title)){
                    // cannot use the following method, spectrum title is empty
                    // !SpectraInput.spectra2targeted_peptide.containsKey(spectrum.getSpectrumTitle())){
                continue;
            }

            if (spectrum.getPrecursor().possibleCharges.length<=0) {
                Cloger.getInstance().logger.warn("There is no precursor charge. The spectrum is ignored!");
            } else {
                charge = spectrum.getPrecursor().possibleCharges[0];
                mass = spectrum.getPrecursor().getMass(charge);
                if( (mass >= massRange[0]) && (mass <= massRange[1]) && (spectrum.getNPeaks() >= minPeaks)){
                    spectrum.spectrumTitle = title;
                    spectrums.add(spectrum);
                }
            }

        }

        if(isZip){
            File DF = new File(ms_file);
            DF.delete();
        }

        return(spectrums);
    }
}
