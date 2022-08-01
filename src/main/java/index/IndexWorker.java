package main.java.index;

import com.compomics.util.experiment.biology.ions.impl.ElementaryIon;
import com.compomics.util.experiment.io.mass_spectrometry.mgf.MgfFileIterator;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import com.compomics.util.gui.waiting.waitinghandlers.WaitingHandlerCLIImpl;
import com.compomics.util.waiting.WaitingHandler;
import main.java.util.Cloger;
import org.tukaani.xz.XZInputStream;
import umich.ms.datatypes.LCMSDataSubset;
import umich.ms.datatypes.lcmsrun.LCMSRunInfo;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scan.StorageStrategy;
import umich.ms.datatypes.scancollection.IScanCollection;
import umich.ms.datatypes.scancollection.impl.ScanCollectionDefault;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.mzml.MZMLFile;
import umich.ms.fileio.filetypes.mzml.MZMLIndex;
import umich.ms.fileio.filetypes.mzxml.MZXMLFile;
import umich.ms.fileio.filetypes.mzxml.MZXMLIndex;
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.zip.GZIPInputStream;

import static main.java.index.BuildMSlibrary.getIntMass;

public final class IndexWorker implements Runnable{


    private String msfile;
    private String outdir;
    private int fold = 10;

    public IndexWorker(String msms_file, String outdir, int fold){
        this.msfile = msms_file;
        this.outdir = outdir;
        this.fold = fold;
        File OD = new File(this.outdir);
        if(!OD.isDirectory()){
            OD.mkdirs();
        }
    }

    @Override
    public void run() {

        long start_time = System.currentTimeMillis();
        Cloger.getInstance().logger.info(Thread.currentThread().getName()+"\t"+this.msfile);

        File F = new File(this.msfile);
        String fullMsFilePath = F.getAbsolutePath();
        String msFileName = F.getName();
        String msFileName_without_suffix = msFileName;
        msFileName_without_suffix = msFileName_without_suffix.replaceAll(".gz$","");
        if(msFileName_without_suffix.toLowerCase().endsWith(".mgf")){
            msFileName_without_suffix = msFileName_without_suffix.replaceAll(".mgf$","");
        }else if(msFileName_without_suffix.toLowerCase().endsWith(".mzml")){
            msFileName_without_suffix = msFileName_without_suffix.replaceAll(".mz[Mm][Ll]$","");
        }

        int cores = 1;

        boolean isGZ = false;
        if (fullMsFilePath.endsWith("gz")) {
            File rawF = new File(fullMsFilePath);
            File tempFile = null;
            try {
                tempFile = File.createTempFile("PepQuery_", rawF.getName().replaceAll(".gz$", ""));
            } catch (IOException e) {
                e.printStackTrace();
            }
            try {
                decompressGzipFile(fullMsFilePath, tempFile.getAbsolutePath());
            } catch (IOException e) {
                e.printStackTrace();
            }
            fullMsFilePath = tempFile.getAbsolutePath();
            isGZ = true;

        }

        int total_spectra = 0;
        int n_msms_without_precursor_charge = 0;

        if(fullMsFilePath.toLowerCase().endsWith(".mzml")) {

            MZMLFile source = null;
            if (fullMsFilePath.endsWith("mzML") || fullMsFilePath.endsWith("mzml")) {
                // System.out.println(path.toString());
                source = new MZMLFile(fullMsFilePath);
            }

            LCMSRunInfo lcmsRunInfo = null;
            try {
                lcmsRunInfo = source.fetchRunInfo();
            } catch (FileParsingException e) {
                e.printStackTrace();
            }
            Cloger.getInstance().logger.info(lcmsRunInfo.toString());
            source.setNumThreadsForParsing(cores);

            MZMLIndex mzMLindex = null;
            try {
                mzMLindex = source.fetchIndex();
            } catch (FileParsingException e) {
                e.printStackTrace();
            }

            if (mzMLindex.size() > 0) {

            } else {
                System.err.println("Parsed index was empty!");
            }

            IScanCollection scans;

            scans = new ScanCollectionDefault(true);
            scans.setDataSource(source);
            try {
                scans.loadData(LCMSDataSubset.MS2_WITH_SPECTRA, StorageStrategy.STRONG);
            } catch (FileParsingException e) {
                e.printStackTrace();
            }


            TreeMap<Integer, IScan> num2scanMap = scans.getMapNum2scan();
            Set<Map.Entry<Integer, IScan>> num2scanEntries = num2scanMap.entrySet();



            for (Map.Entry<Integer, IScan> next : num2scanEntries) {
                IScan scan = next.getValue();
                if (scan.getSpectrum() != null) {
                    //System.out.println(scan.getNum());
                    if (scan.getMsLevel() == 2) {
                        total_spectra++;
                        // System.out.println(scan.getNum()+":"+scan.getPrecursor().getCharge()+":"+scan.getPrecursor().getMzTarget());
                        if(scan.getPrecursor().getCharge() != null){
                            int charge = scan.getPrecursor().getCharge();
                            // scan number is original scan number
                            //double mz = scan.getPrecursor().getMzTarget();
                            double mz = scan.getPrecursor().getMzTargetMono();
                            double mass = getMass(mz, charge);
                            try {
                                save2file(mass, this.fold, asMgf(scan, msFileName_without_suffix), outdir);
                            } catch (IOException e) {
                                e.printStackTrace();
                            }
                        }else{
                            n_msms_without_precursor_charge++;
                            Cloger.getInstance().logger.warn("Precursor charge is empty: file=>" + msFileName_without_suffix + ", scan="+scan.getNum());
                            // if charge is invalid
                            int charge = 2;
                            // scan number is original scan number
                            // double mz = scan.getPrecursor().getMzTarget();
                            double mz = scan.getPrecursor().getMzTargetMono();
                            double mass = getMass(mz, charge);
                            try {
                                save2file(mass, this.fold, asMgf(scan, msFileName_without_suffix, charge), outdir);
                            } catch (IOException e) {
                                e.printStackTrace();
                            }

                            charge = 3;
                            // scan number is original scan number
                            // mz = scan.getPrecursor().getMzTarget();
                            mz = scan.getPrecursor().getMzTargetMono();
                            mass = getMass(mz, charge);
                            try {
                                save2file(mass, this.fold, asMgf(scan, msFileName_without_suffix, charge), outdir);
                            } catch (IOException e) {
                                e.printStackTrace();
                            }

                        }
                    }
                }
            }
        }else if(fullMsFilePath.toLowerCase().endsWith(".mzxml")){
            MZXMLFile source = null;
            if (fullMsFilePath.endsWith("mzXML") || fullMsFilePath.endsWith("mzxml")) {
                // System.out.println(path.toString());
                source = new MZXMLFile(fullMsFilePath);
            }

            LCMSRunInfo lcmsRunInfo = null;
            try {
                lcmsRunInfo = source.fetchRunInfo();
            } catch (FileParsingException e) {
                e.printStackTrace();
            }
            Cloger.getInstance().logger.info(lcmsRunInfo.toString());
            source.setNumThreadsForParsing(cores);

            MZXMLIndex mzxmlIndex = null;
            try {
                mzxmlIndex = source.fetchIndex();
            } catch (FileParsingException e) {
                e.printStackTrace();
            }

            if (mzxmlIndex.size() > 0) {

            } else {
                System.err.println("Parsed index was empty!");
            }

            IScanCollection scans;

            scans = new ScanCollectionDefault(true);
            scans.setDataSource(source);
            try {
                scans.loadData(LCMSDataSubset.MS2_WITH_SPECTRA, StorageStrategy.STRONG);
            } catch (FileParsingException e) {
                e.printStackTrace();
            }


            TreeMap<Integer, IScan> num2scanMap = scans.getMapNum2scan();
            Set<Map.Entry<Integer, IScan>> num2scanEntries = num2scanMap.entrySet();



            for (Map.Entry<Integer, IScan> next : num2scanEntries) {
                IScan scan = next.getValue();
                if (scan.getSpectrum() != null) {
                    //System.out.println(scan.getNum());
                    if (scan.getMsLevel() == 2) {
                        total_spectra++;
                        // System.out.println(scan.getNum()+":"+scan.getPrecursor().getCharge()+":"+scan.getPrecursor().getMzTarget());
                        if(scan.getPrecursor().getCharge() != null){
                            int charge = scan.getPrecursor().getCharge();
                            // scan number is original scan number
                            //double mz = scan.getPrecursor().getMzTarget();
                            // for mzXML, there is no getMzTargetMono
                            double mz = scan.getPrecursor().getMzTarget();
                            double mass = getMass(mz, charge);
                            try {
                                save2file(mass, this.fold, asMgf(scan, msFileName_without_suffix), outdir);
                            } catch (IOException e) {
                                e.printStackTrace();
                            }
                        }else{
                            n_msms_without_precursor_charge++;
                            Cloger.getInstance().logger.warn("Precursor charge is empty: file=>" + msFileName_without_suffix + ", scan="+scan.getNum());
                            // if charge is invalid
                            int charge = 2;
                            // scan number is original scan number
                            // double mz = scan.getPrecursor().getMzTarget();
                            double mz = scan.getPrecursor().getMzTargetMono();
                            double mass = getMass(mz, charge);
                            try {
                                save2file(mass, this.fold, asMgf(scan, msFileName_without_suffix, charge), outdir);
                            } catch (IOException e) {
                                e.printStackTrace();
                            }

                            charge = 3;
                            // scan number is original scan number
                            // mz = scan.getPrecursor().getMzTarget();
                            mz = scan.getPrecursor().getMzTargetMono();
                            mass = getMass(mz, charge);
                            try {
                                save2file(mass, this.fold, asMgf(scan, msFileName_without_suffix, charge), outdir);
                            } catch (IOException e) {
                                e.printStackTrace();
                            }

                        }
                    }
                }
            }
        }else if(fullMsFilePath.toLowerCase().endsWith(".mgf")){

            File mgfFile = new File(fullMsFilePath);
            WaitingHandler waitingHandler = new WaitingHandlerCLIImpl();
            waitingHandler.setDisplayProgress(false);

            MgfFileIterator mgfFileIterator = new MgfFileIterator(mgfFile, waitingHandler);
            String title;
            Spectrum spectrum = new Spectrum();
            double mass;
            int charge;
            long spectrum_index = 0;
            String spectrumTitle;
            while ((title = mgfFileIterator.next()) != null) {

                spectrum = mgfFileIterator.getSpectrum();

                spectrum_index++;
                total_spectra++;

                // scan number is not original scan number any more
                String scanNumber = String.valueOf(spectrum_index);

                if (spectrum.getPrecursor().possibleCharges.length<=0) {
                    n_msms_without_precursor_charge++;

                    Cloger.getInstance().logger.warn("Precursor charge is empty: " + spectrum.getSpectrumTitle());
                    // must process spectra without precursor charge information
                    // will consider 2+ and 3+ in default
                    charge = 2;
                    mass = spectrum.getPrecursor().getMass(charge);
                    spectrumTitle = msFileName_without_suffix+":"+scanNumber+":"+charge;
                    try {
                        save2file(mass,this.fold,asMgf(spectrum, spectrumTitle,charge,scanNumber),outdir);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }

                    charge = 3;
                    mass = spectrum.getPrecursor().getMass(charge);
                    try {
                        save2file(mass,this.fold,asMgf(spectrum, spectrumTitle,charge,scanNumber),outdir);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }

                } else {

                    charge = spectrum.getPrecursor().possibleCharges[0];
                    mass = spectrum.getPrecursor().getMass(charge);

                    // Need to process spectra without precursor charge.

                    try {
                        save2file(mass,this.fold,asMgf(spectrum, msFileName_without_suffix,scanNumber),outdir);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }

            }
        }

        if (isGZ) {
            Cloger.getInstance().logger.warn("Delete " + fullMsFilePath);
            try {
                Files.delete(Paths.get(fullMsFilePath));
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        Cloger.getInstance().logger.info(msFileName + "\t" + total_spectra+"\tspectra without precursor charge:"+n_msms_without_precursor_charge+"("+String.format("%.4f%%)",100.0*n_msms_without_precursor_charge/total_spectra));
        long ctime = System.currentTimeMillis();
        double t = 1.0*(ctime  - start_time)/1000.0/60.0;
        String out = String.format("%.2f",t);
        out = out + " min";

        try {
            BufferedWriter nWriter = new BufferedWriter(new FileWriter(new File(this.outdir+"/count.txt")));
            nWriter.write(total_spectra +"\n");
            nWriter.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        Cloger.getInstance().logger.info(Thread.currentThread().getName()+"\t"+this.msfile+"\t"+this.outdir+"\t"+out);

    }

    public static void decompressGzipFile(String input, String output) throws IOException {

        InputStream fileStream = new FileInputStream(input);
        InputStream gzipStream = new GZIPInputStream(fileStream);
        Reader decoder = new InputStreamReader(gzipStream);
        BufferedReader gReader = new BufferedReader(decoder);

        BufferedWriter pWriter = new BufferedWriter(new FileWriter(new File(output)));
        String line;
        while((line = gReader.readLine())!=null){
            pWriter.write(line+"\n");
        }
        gReader.close();
        pWriter.close();

    }

    public static void decompressXzipFile(String input, String output) throws IOException {

        InputStream fileStream = new FileInputStream(input);
        InputStream gzipStream = new XZInputStream(fileStream);
        Reader decoder = new InputStreamReader(gzipStream);
        BufferedReader gReader = new BufferedReader(decoder);

        BufferedWriter pWriter = new BufferedWriter(new FileWriter(new File(output)));
        String line;
        while((line = gReader.readLine())!=null){
            pWriter.write(line+"\n");
        }
        gReader.close();
        pWriter.close();

    }

    private String asMgf(IScan scan, String filename){

        int charge = scan.getPrecursor().getCharge();
        String out = asMgf(scan,filename,charge);
        return out;
    }

    private String asMgf(IScan scan, String filename, int charge){

        // minute
        double rt = scan.getRt();
        double intensity = 0;
        if(scan.getPrecursor().getIntensity() != null){
            intensity = scan.getPrecursor().getIntensity();
        }
        //double mz = scan.getPrecursor().getMzTarget();
        double mz;
        if(scan.getPrecursor().getMzTargetMono()!=null){
            mz = scan.getPrecursor().getMzTargetMono();
        }else{
            mz = scan.getPrecursor().getMzTarget();
        }

        // System.out.println(scan.getPrecursor().getMzTarget()+"\t"+scan.getPrecursor().getMzTargetMono());
        double[] intArray = scan.getSpectrum().getIntensities();
        double[] mzArray = scan.getSpectrum().getMZs();
        int scan_number = scan.getNum();

        StringBuilder stringBuilder = new StringBuilder();

        stringBuilder.append("BEGIN IONS\n");
        stringBuilder.append("TITLE=").append(filename).append(":").append(scan_number).append(":").append(charge).append("\n");
        stringBuilder.append("PEPMASS=").append(mz).append(" ").append(intensity).append("\n");
        // is there a way to check rt unit from scan.getRt()?
        stringBuilder.append("RTINSECONDS=").append(rt*60).append("\n");
        stringBuilder.append("CHARGE=").append(charge).append("+\n");
        stringBuilder.append("SCANS=").append(scan_number).append("\n");
        for(int i=0;i<intArray.length;i++){
            stringBuilder.append(String.format("%.4f",mzArray[i])).append(" ").append(String.format("%.2f",intArray[i])).append("\n");
        }
        stringBuilder.append("END IONS\n");
        return(stringBuilder.toString());

    }

    public static String asMgf(Spectrum spectrum, String filename, String scan_number){

        int charge = spectrum.getPrecursor().possibleCharges[0];
        String spectrumTitle = filename+":"+scan_number+":"+charge;
        String res = asMgf(spectrum,spectrumTitle,charge, scan_number);
        return(res);

    }

    public static  String asMgf(Spectrum spectrum, String spectrumTitle,int charge, String scan_number){

        double intensity = spectrum.getPrecursor().intensity;
        double mz = spectrum.getPrecursor().mz;
        //HashMap<Double, Peak> peakMap = spectrum.getPeakMap();
        double[] mzArray = spectrum.mz;
        double[] intensityArray = spectrum.intensity;

        StringBuilder stringBuilder = new StringBuilder();

        stringBuilder.append("BEGIN IONS\n");
        stringBuilder.append("TITLE=").append(spectrumTitle).append("\n");
        stringBuilder.append("PEPMASS=").append(mz).append(" ").append(intensity).append("\n");
        if(spectrum.getPrecursor().rt >= 0.0) {
            double rt = spectrum.getPrecursor().rt;
            stringBuilder.append("RTINSECONDS=").append(rt).append("\n");
        }
        stringBuilder.append("CHARGE=").append(charge).append("+\n");
        if(scan_number != null){
            stringBuilder.append("SCANS=").append(scan_number).append("\n");
        }
        for(int i=0;i<mzArray.length;i++){
            stringBuilder.append(String.format("%.4f",mzArray[i])).append(" ").append(String.format("%.2f",intensityArray[i])).append("\n");
        }
        stringBuilder.append("END IONS\n");
        return(stringBuilder.toString());

    }



    /**
     * Returns the mass of the precursor with the given charge.
     *
     * @param chargeValue the value of the charge
     *
     * @return the mass of the precursor with the given charge
     */
    private double getMass(double mz, int chargeValue) {
        return mz * chargeValue - chargeValue * ElementaryIon.proton.getTheoreticMass();
    }



    private void save2file(double mass, int fold, String spectrum, String outdir) throws IOException {

        int intMass = getIntMass(mass,fold);
        String out_file = outdir + File.separator + intMass + ".mgf";
        File F = new File(out_file);
        if(!F.isFile()){
            F.createNewFile();
        }
        FileWriter fw = new FileWriter(F,true);
        BufferedWriter bw = new BufferedWriter(fw);
        bw.write(spectrum+"\n");
        bw.close();
    }
}
