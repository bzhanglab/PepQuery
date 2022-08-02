package main.java.msio;

import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import main.java.pg.CParameter;
import main.java.pg.PeptideInput;
import main.java.util.Cloger;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import static main.java.pg.SpectraInput.getRangeOfMass;

public class MsdataMatch {




    public MsdataMatch(){

    }

    public HashSet<String> get_all_indexed_ms_file(String msms_library_path, ArrayList<PeptideInput> peptideInputs, String out_dir, int n_cpu){

        // The process to get candidate index files should be the same with that used in MsLibrarySearchWorker

        msms_library_path = msms_library_path.replaceFirst("/$","");

        HashSet<String> indexed_ms_files = new HashSet<>();

        for(PeptideInput peptideInput: peptideInputs) {
            // key: peptide mass, value: peptide form in peptideInput.getPtmIsoforms()
            // This is used to speed up spectra finding. For the peptide forms with same mass, only one time finding is
            // needed.
            HashMap<Double, ArrayList<Integer>> pepMass2index = new HashMap<>();
            for (int j = 0; j < peptideInput.getPtmIsoforms().size(); j++) {
                // for each peptide form, get spectra
                ArrayList<Spectrum> rsSpectra = new ArrayList<>();
                double pepIsoformMass = peptideInput.getPtmIsoforms().get(j).getMass();
                if (pepMass2index.containsKey(pepIsoformMass)) {
                    pepMass2index.get(pepIsoformMass).add(j);
                } else {
                    ArrayList<Integer> indexList = new ArrayList<>();
                    indexList.add(j);
                    pepMass2index.put(pepIsoformMass, indexList);
                }
            }

            //for (int j = 0; j < peptideInput.getPtmIsoforms().size(); j++) {
            for (double pepIsoformMass : pepMass2index.keySet()) {
                // for each unique peptide mass, get spectra
                ArrayList<Spectrum> rsSpectra = new ArrayList<>();
                // double pepIsoformMass = peptideInput.getPtmIsoforms().get(j).getMass();
                ArrayList<Double> isotope_masses = CParameter.get_peptide_mass_with_isotope_error(pepIsoformMass);
                for (double i_mass : isotope_masses) {


                    double[] massRange = getRangeOfMass(i_mass, CParameter.tol, CParameter.tolu.equalsIgnoreCase("ppm"));
                    int left_range = (int) Math.round(massRange[0] * 10);
                    int right_range = (int) Math.round(massRange[1] * 10);

                    //System.out.println(i_mass+"\t"+massRange[0]+"\t"+massRange[1]);

                    for (int i = left_range; i <= right_range; i++) {
                        String i_mgf;
                        if(msms_library_path.toLowerCase().startsWith("s3:")) {
                            i_mgf = msms_library_path + "/" + i + ".mgf";
                        }else{
                            i_mgf = msms_library_path + File.separator + i + ".mgf";
                        }
                        indexed_ms_files.add(i_mgf);
                    }
                }
            }
        }

        File OD = new File(out_dir);
        if(!OD.isDirectory()){
            OD.mkdirs();
        }
        // download all required index files
        if(n_cpu > indexed_ms_files.size() && indexed_ms_files.size() > 0){
            n_cpu = indexed_ms_files.size();
        }
        Cloger.getInstance().logger.info("Used CPUs: "+n_cpu);
        Cloger.getInstance().logger.info("Download "+ indexed_ms_files.size() + " index MS/MS files ...");

        ExecutorService fixedThreadPoolScore = Executors.newFixedThreadPool(n_cpu);

        // For each target peptide, retrieve its matched spectra from the spectra library.
        for (String file: indexed_ms_files) {
            fixedThreadPoolScore.execute(new IndexFileDownloadWorker(file,out_dir));

        }

        fixedThreadPoolScore.shutdown();
        try {
            fixedThreadPoolScore.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        }

        return indexed_ms_files;

    }
}
