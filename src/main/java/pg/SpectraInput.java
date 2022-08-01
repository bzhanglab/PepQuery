package main.java.pg;

import com.compomics.util.experiment.io.mass_spectrometry.mgf.MgfFileIterator;
import com.compomics.util.experiment.mass_spectrometry.spectra.Precursor;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;

import com.compomics.util.gui.waiting.waitinghandlers.WaitingHandlerCLIImpl;
import com.compomics.util.waiting.WaitingHandler;
import main.java.PSMMatch.JPeptide;
import main.java.msio.MsLibrarySearchWorker;
import main.java.msio.MsdataMatch;
import main.java.util.Cloger;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.sql.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Read spectrum and filter according to peptide sequence
 */
public class SpectraInput {

    /**
     * It is used to store index file downloaded from s3.
     */
    public static String out_dir = "./";

    public static boolean isTargetPeptide2spectrumScored = false;

    /**
     * All matched spectra to any target peptide forms.
     */
    public static HashMap<String, Spectrum> spectraMap = new HashMap<>();

    /**
     * Peptide spectrum match from input file (-pep). If this is not empty, PepQuery will only compare peptides against
     * the spectra from the input file (-pep) and ignore all other spectra from MS/MS data file provided.
     * This is useful when users are only interested in validating specific PSMs matches from other peptide
     * identification method. It is also used to record the status of each targeted PSM provided by users.
     */
    public static HashMap<String, HashMap<String,TargetedPSM>> pep2targeted_spectra = new HashMap<>();

    /**
     * Peptide spectrum match from input file (-pep). If this is not empty, PepQuery will only compare peptides against
     * the spectra from the input file (-pep) and ignore all other spectra from MS/MS data file provided.
     * This is useful when users are only interested in validating specific PSMs matches from other peptide
     * identification method.
     */
    public static HashMap<String, HashMap<String, TargetedPSM>> spectra2targeted_peptide = new HashMap<>();

    /**
     * Remove all spectra
     */
    public static void clear(){
        spectraMap.clear();
        pep2targeted_spectra.clear();
        spectra2targeted_peptide.clear();
    }



    // used for testing
    public static void main(String[] args) {


    }


    public static Spectrum getSpectrum(String spectrum_title){
        return(spectraMap.get(spectrum_title));
    }

    // When taking multiple peptides as input, build a peptide index for fast searching

    /**
     *
     * @param spectraFile
     * @param peptideInputs
     * @throws IOException
     * @throws InterruptedException
     */
    public static void readMSMSfast(String spectraFile, ArrayList<PeptideInput> peptideInputs) throws IOException, InterruptedException {


        // bin=0.1Da

        ConcurrentHashMap<Integer,ArrayList<JPeptide>> pIndex = new ConcurrentHashMap<>();
        for(PeptideInput peptideInput: peptideInputs){
            for(JPeptide jPeptide2: peptideInput.getPtmIsoforms()){
                //rmass = (int) Math.round(10.0*jPeptide2.peptide.getMass());
                int rmass = (int) Math.round(10.0*jPeptide2.getMass());
                if(pIndex.containsKey(rmass)){
                    pIndex.get(rmass).add(jPeptide2);
                }else{
                    ArrayList<JPeptide> jps = new ArrayList<>();
                    jps.add(jPeptide2);
                    pIndex.put(rmass,jps);
                }
            }
        }

        Cloger.getInstance().logger.info("Build target peptide index done.");

        boolean isS3 = false;

        if(spectraFile.startsWith("s3:")){
            isS3 = true;
            String objPath = spectraFile;
            File TF = new File(objPath);
            String file_name = TF.getName();
            if(file_name.endsWith(".xz")){
                file_name = file_name.replaceFirst(".xz","");
            }else if(file_name.endsWith(".gz")){
                file_name = file_name.replaceFirst(".gz","");
            }
            File tempFile = null;
            try {
                tempFile = File.createTempFile(Thread.currentThread().getName(), "_"+file_name);
            } catch (IOException e) {
                e.printStackTrace();
            }

            spectraFile = tempFile.getAbsolutePath();

            if(S3Interface.getInstance().download(objPath,spectraFile)){
                Cloger.getInstance().logger.info("Use MS/MS data from S3: "+objPath);
            }else{
                Cloger.getInstance().logger.info(objPath+ " doesn't exist");
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ExecutorService fixedThreadPoolScore = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());

        //SpectrumFactory spectrumFactory = SpectrumFactory.getInstance();
        File mgfFile = new File(spectraFile);
        //spectrumFactory.addSpectra(mgfFile, null);
        WaitingHandler waitingHandler = new WaitingHandlerCLIImpl();
        waitingHandler.setDisplayProgress(false);

        MgfFileIterator mgfFileIterator = new MgfFileIterator(mgfFile, waitingHandler);
        String title;
        Spectrum spectrum;
        int i = 0;
        while ((title = mgfFileIterator.next()) != null) {

            i++;
            spectrum = mgfFileIterator.getSpectrum();

            //System.out.println(title);
            //System.out.println(CParameter.indexType);
            if(CParameter.indexType != 2) {
                spectrum.spectrumTitle = String.valueOf(i);
                //System.out.println(CParameter.tagIndexType+"\t"+spectrum.getSpectrumTitle());
            }else{
                spectrum.spectrumTitle = title;
            }

            // This is for the mgf file from indexed MS/MS spectra library, so each spectrum in the mgf file should have
            // a valid and unique spectrum title.

            // If SpectraInput.spectra2targeted_peptide is not empty, it means the input for -pep contains spectra for
            // each targeted peptide. If this is the case, we only need to score peptides against the targeted spectra
            if(SpectraInput.spectra2targeted_peptide.size()>=1 &&
                    !SpectraInput.spectra2targeted_peptide.containsKey(spectrum.getSpectrumTitle())){
                continue;
            }

            //else{
                //System.out.println(spectrum.getSpectrumTitle());
            //}
            //MSnSpectrum spectrum = (MSnSpectrum) spectrumFactory.getSpectrum(msfileName, String.valueOf(i));
            // if the peak number is less than minPeaks, then this spectra will be removed
            if(spectrum.getNPeaks()<=CParameter.minPeaks){
                continue;
            }
            // 2018-04-19
            // It's possible that the spectrum doesn't have charge information.
            fixedThreadPoolScore.execute(new FindCandidateSpectraAndScoreWorker(pIndex,spectrum));

            if( (i%10000)==0){
                Cloger.getInstance().logger.info("Finished:"+i);

            }

        }

        fixedThreadPoolScore.shutdown();
        fixedThreadPoolScore.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);


        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        int npsms = 0;

        // Some spectra may not have charge state information. Save the charge state information used for scoring
        // for these spectra.
        HashMap<String,Integer> spectrumID2charge = new HashMap<>();
        HashSet<String> matchSpectra = new HashSet<>(128);
        for(PeptideInput peptideInput: peptideInputs){
            for(JPeptide jPeptide2: peptideInput.getPtmIsoforms()){
                if(!jPeptide2.spectraIndexs.isEmpty()){
                    matchSpectra.addAll(jPeptide2.spectraIndexs);
                    npsms = npsms + jPeptide2.spectraIndexs.size();
                    for(String spectrumID : jPeptide2.spectrumID2charge.keySet()){
                        spectrumID2charge.put(spectrumID,jPeptide2.spectrumID2charge.get(spectrumID));
                    }
                }
            }
        }

        Cloger.getInstance().logger.info("Matched spectra: "+matchSpectra.size());
        Cloger.getInstance().logger.info("Matched PSMs: "+npsms);
        if(spectrumID2charge.size()>=1) {
            Cloger.getInstance().logger.info("Spectra without charge information: " + spectrumID2charge.size());
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        waitingHandler = new WaitingHandlerCLIImpl();
        waitingHandler.setDisplayProgress(false);

        mgfFileIterator = new MgfFileIterator(mgfFile, waitingHandler);
        i=0;
        while ((title = mgfFileIterator.next()) != null) {

            i++;
            spectrum = mgfFileIterator.getSpectrum();

            if(CParameter.indexType != 2) {
                spectrum.spectrumTitle = String.valueOf(i);
            }else{
                spectrum.spectrumTitle = title;
            }
            //MSnSpectrum spectrum = (MSnSpectrum) spectrumFactory.getSpectrum(msfileName, String.valueOf(i));
            // if the peak number is less than minPeaks, then this spectra will be removed
            if(matchSpectra.contains(spectrum.getSpectrumTitle())){
                if(spectrum.getPrecursor().possibleCharges.length<=0){
                    // no charge information
                    // We need to change the precursor of this spectrum
                    //ArrayList<Integer> possibleCharges = new ArrayList<>();
                    int charge = spectrumID2charge.get(spectrum.getSpectrumTitle());
                    int possibleCharges [] = new int[]{charge};
                    Precursor precursor = new Precursor(0.0,spectrum.getPrecursor().mz,possibleCharges);
                    spectrum.precursor = precursor;
                }
                spectraMap.put(spectrum.getSpectrumTitle(),spectrum);
            }

            if( (i%10000)==0){
                Cloger.getInstance().logger.info("Finished:"+i);

            }

        }

        pIndex = null;
        mgfFileIterator = null;
        matchSpectra = null;

        if(isS3){
            File FF = new File(spectraFile);
            FF.delete();
        }

        Cloger.getInstance().logger.info("Saved spectra: "+spectraMap.size());

        Cloger.getInstance().logger.info("Get matched spectra done!");

        isTargetPeptide2spectrumScored = true;


    }


    public static void readMSMSdb2(String sqldb, ArrayList<PeptideInput> peptideInputs) throws IOException, SQLException, ClassNotFoundException, InterruptedException {

        Connection connection = DriverManager.getConnection("jdbc:sqlite:" + sqldb);
        PreparedStatement pstmt = connection.prepareStatement("select id,spectrum from msmsdb where mass >= ? and mass <= ?");

        int matchSpectra = 0;

        PeptideInput peptideInput;
        double peptideMass;
        double pepIsoformMass;

        // don't use multiple threads here. Multiple threads will cause problem here.
        ExecutorService fixedThreadPoolScore = Executors.newFixedThreadPool(1);

        ConcurrentHashMap<String,Spectrum> matchedSpectra = new ConcurrentHashMap<>();

        int psize = peptideInputs.size();
        for (int pindex=0;pindex < psize;pindex++) {
            peptideInput = peptideInputs.get(pindex);
            for (int j = 0; j < peptideInput.getPtmIsoforms().size(); j++) {
                pepIsoformMass = peptideInput.getPtmIsoforms().get(j).getMass();
                double massRange2[] = getRangeOfMass(pepIsoformMass,CParameter.tol,CParameter.tolu.equalsIgnoreCase("ppm"));
                pstmt.setDouble(1,massRange2[0]);
                pstmt.setDouble(2,massRange2[1]);

                ResultSet rs = pstmt.executeQuery();
                ArrayList<Spectrum> rsSpectra2 = getSpectraFromSQL(rs,CParameter.minPeaks);
                for(Spectrum mSnSpectrum: rsSpectra2){
                    fixedThreadPoolScore.execute(new Peptide2SpectrumScoreWorker(peptideInput.getPtmIsoforms().get(j),mSnSpectrum, matchedSpectra));

                }

                matchSpectra = matchSpectra + rsSpectra2.size();
            }
        }

        fixedThreadPoolScore.shutdown();
        fixedThreadPoolScore.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);



        for(String title:matchedSpectra.keySet()){
            spectraMap.put(title,matchedSpectra.get(title));
        }

        //matchedSpectra.clear();
        Cloger.getInstance().logger.info("Matched spectra: "+spectraMap.size());
        Cloger.getInstance().logger.info("Matched PSMs: "+matchSpectra);

        pstmt.close();
        connection.close();

        matchedSpectra = null;

        isTargetPeptide2spectrumScored = true;

    }

    /**
     * Retrieve spectra for target peptides from indexed MS/MS data
     * @param msms_library_path Indexed MS/MS data
     * @param peptideInputs Target peptide sequences
     * @throws InterruptedException
     */
    public static void readSpectraFromMSMSlibrary(String msms_library_path, ArrayList<PeptideInput> peptideInputs) throws InterruptedException {

        if(msms_library_path.startsWith("s3:")){
            S3Interface.getInstance();
        }

        int ncpu = CParameter.cpu;
        if( ncpu == 0){
            ncpu = Runtime.getRuntime().availableProcessors();
        }

        // The MS/MS index files downloaded to local computer
        HashSet<String> localTempIndexFiles = new HashSet<>();
        if(msms_library_path.toLowerCase().startsWith("s3://")) {
            MsdataMatch msdataMatch = new MsdataMatch();
            String local_index_file_dir = out_dir + File.separator + "index";
            localTempIndexFiles = msdataMatch.get_all_indexed_ms_file(msms_library_path, peptideInputs, local_index_file_dir, ncpu);
            msms_library_path = local_index_file_dir;
        }

        int psize = peptideInputs.size();

        if(ncpu > psize && psize > 0){
            ncpu = psize;
        }
        Cloger.getInstance().logger.info("Used CPUs: "+ncpu);

        ExecutorService fixedThreadPoolScore = Executors.newFixedThreadPool(ncpu);

        // For each target peptide, retrieve its matched spectra from the spectra library.
        for (int pindex=0;pindex < psize;pindex++) {
            PeptideInput peptideInput = peptideInputs.get(pindex);
            // Find candidate spectra for each target peptide
            fixedThreadPoolScore.execute(new MsLibrarySearchWorker(peptideInput,msms_library_path));

        }

        fixedThreadPoolScore.shutdown();
        fixedThreadPoolScore.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);

        // save all matched spectra to spectraMap: spectrum title -> Spectrum object
        for (int pindex=0;pindex < psize;pindex++) {
            spectraMap.putAll(peptideInputs.get(pindex).getSpectra());
            peptideInputs.get(pindex).clearSpectra();

        }

        Cloger.getInstance().logger.info("Matched spectra: "+spectraMap.size());
        //System.out.println("Matched PSMs: "+matchSpectra);

        isTargetPeptide2spectrumScored = true;

        if(msms_library_path.startsWith("s3:")){
            //S3Interface.getInstance().shutdown();

        }

        // The MS/MS index files downloaded to local computer
        if(localTempIndexFiles.size()>0){
            // delete downloaded index files to save storage
            Cloger.getInstance().logger.info("Delete downloaded MS/MS index files.");
            for(String file: localTempIndexFiles){
                File F = new File(file);
                File localF = new File(msms_library_path+File.separator+F.getName());
                if(localF.isFile()){
                    localF.delete();
                }
            }
            // delete the folder
            File IndexF = new File(msms_library_path);
            if(IndexF.isDirectory()){
                IndexF.delete();
            }
        }

    }


    /**
     * Calculate the left and right values for specified tol.
     * @param pepmass peptide mass
     * @param tol   tol
     * @param isPpm true if the unit of tol is ppm.
     * @return An array with two elements.
     */
    public static double[] getRangeOfMass(double pepmass, double tol, boolean isPpm){
        double [] massRange = new double [2];
        if(isPpm){
            double x = 1.0*(tol * pepmass) / (1.0e6);
            massRange[0] = pepmass-x;
            massRange[1] = pepmass+x;

        }else{
            massRange[0] = pepmass-tol;
            massRange[1] = pepmass+tol;
        }
        return massRange;

    }



    public static ArrayList<Spectrum> getSpectraFromSQL(ResultSet rs, int minPeaks) throws SQLException, IOException, ClassNotFoundException {
        ArrayList<Spectrum> spectrums = new ArrayList<>();

        while (rs.next()) {

            byte buf[] = rs.getBytes("spectrum");
            String id = rs.getString("id");
            ObjectInputStream objectInputStream = new ObjectInputStream(new ByteArrayInputStream(buf));
            Spectrum p = (Spectrum) objectInputStream.readObject();
            if(p.getNPeaks() > minPeaks){
                p.spectrumTitle = id;
                spectrums.add(p);
            }

        }
        return spectrums;
    }


}
