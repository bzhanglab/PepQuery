package main.java.plot;

import com.compomics.util.experiment.biology.ions.Ion;
import com.compomics.util.experiment.biology.ions.NeutralLoss;
import com.compomics.util.experiment.biology.ions.impl.PeptideFragmentIon;
import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.identification.matches.IonMatch;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.experiment.identification.spectrum_annotation.AnnotationParameters;
import com.compomics.util.experiment.identification.spectrum_annotation.SpecificAnnotationParameters;
import com.compomics.util.experiment.identification.spectrum_annotation.SpectrumAnnotator;
import com.compomics.util.experiment.identification.spectrum_annotation.spectrum_annotators.PeptideSpectrumAnnotator;
import com.compomics.util.experiment.identification.spectrum_assumptions.PeptideAssumption;
import com.compomics.util.experiment.io.mass_spectrometry.mgf.MgfFileIterator;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import com.compomics.util.gui.waiting.waitinghandlers.WaitingHandlerCLIImpl;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import com.compomics.util.parameters.identification.search.ModificationParameters;
import com.compomics.util.waiting.WaitingHandler;
import main.java.OpenModificationSearch.ModificationDB;
import main.java.PSMMatch.HyperscoreMatch;
import main.java.PSMMatch.JPeptide;
import main.java.pg.CParameter;
import main.java.pg.JSequenceProvider;
import main.java.pg.PeptideSearchMT;
import main.java.util.Cloger;
import org.apache.commons.cli.*;
import org.apache.commons.lang3.StringUtils;
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

/**
 * Peptide spectrum annotation
 */
public class PeakAnnotation {


    public static double intensityLimit = 0.03;


    public static void main(String[] args) throws IOException, ParseException {

        Cloger.getInstance().logger.info("Start analysis");
        Cloger.getInstance().logger.info(StringUtils.join(args," "));

        HyperscoreMatch.generateFactorialValues(60);

        Options options = new Options();


        options.addOption("ms", true, "Spectrum file used for identification, mgf format");

        options.addOption("i", true, "Peptide identification txt format file");
        options.addOption("o", true, "Output annotation file");
         options.addOption("itol", true, "Fragment ion m/z tolerance, default is 0.6da");

        options.addOption("h", false, "Help");


        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = parser.parse(options, args);
        if (cmd.hasOption("h") || cmd.hasOption("help") || args.length == 0) {
            HelpFormatter f = new HelpFormatter();

            f.printHelp("Options", options);
            System.exit(0);
        }

        String psm_file = cmd.getOptionValue("i");
        String mgf = cmd.getOptionValue("ms");
        String outfile = cmd.getOptionValue("o");
        double itol = 0.6;
        if(cmd.hasOption("itol")){
            itol = Double.valueOf(cmd.getOptionValue("itol"));
        }
        CParameter.itol = itol;
        plotPSM(psm_file,mgf,outfile);
    }

    public static String getPeakAnnotation(Peptide objPeptide, Spectrum spectrum, boolean lossWaterNH3){

        PeptideSpectrumAnnotator peptideSpectrumAnnotator = new PeptideSpectrumAnnotator();

        //String spectrumKey = spectrum.getSpectrumKey();
        int charge = spectrum.getPrecursor().possibleCharges[0];
        PeptideAssumption peptideAssumption = new PeptideAssumption(objPeptide, charge);
        SpecificAnnotationParameters specificAnnotationPreferences = new SpecificAnnotationParameters();

        HashSet<Integer> charges = new HashSet<>(4);
        int precursorCharge = peptideAssumption.getIdentificationCharge();
        if (precursorCharge == 1) {
            charges.add(precursorCharge);
        } else {
            //for (int c = 1; c < precursorCharge; c++) {
            for (int c = 1; c < precursorCharge; c++) {
                charges.add(c);
            }
        }
        specificAnnotationPreferences.setSelectedCharges(charges);

        specificAnnotationPreferences.addIonType(Ion.IonType.PEPTIDE_FRAGMENT_ION, PeptideFragmentIon.B_ION);
        specificAnnotationPreferences.addIonType(Ion.IonType.PEPTIDE_FRAGMENT_ION, PeptideFragmentIon.Y_ION);
        specificAnnotationPreferences.addIonType(Ion.IonType.PRECURSOR_ION);
        specificAnnotationPreferences.addIonType(Ion.IonType.IMMONIUM_ION);

        specificAnnotationPreferences.setFragmentIonAccuracy(CParameter.itol);
        specificAnnotationPreferences.setFragmentIonPpm(false);
        //specificAnnotationPreferences private boolean neutralLossesAuto = true;
        specificAnnotationPreferences.setNeutralLossesAuto(false);
        specificAnnotationPreferences.clearNeutralLosses();
        // this is important
        specificAnnotationPreferences.setPrecursorCharge(precursorCharge);

        if(lossWaterNH3) {
            specificAnnotationPreferences.addNeutralLoss(NeutralLoss.H2O);
            specificAnnotationPreferences.addNeutralLoss(NeutralLoss.NH3);
        }

        AnnotationParameters annotationSettings = new AnnotationParameters();
        annotationSettings.setTiesResolution(SpectrumAnnotator.TiesResolution.mostIntense);
        annotationSettings.setFragmentIonAccuracy(CParameter.itol);
        annotationSettings.setFragmentIonPpm(false);
        annotationSettings.setIntensityLimit(intensityLimit);
        annotationSettings.setNeutralLossesSequenceAuto(false);

        if(lossWaterNH3) {
            annotationSettings.addNeutralLoss(NeutralLoss.H2O);
            annotationSettings.addNeutralLoss(NeutralLoss.NH3);
        }

        if(ModificationDB.getInstance().getModificationString(objPeptide).toLowerCase().contains("phosphorylation")){
            annotationSettings.addNeutralLoss(NeutralLoss.H3PO4);
            annotationSettings.addNeutralLoss(NeutralLoss.HPO3);
            specificAnnotationPreferences.addNeutralLoss(NeutralLoss.H3PO4);
            specificAnnotationPreferences.addNeutralLoss(NeutralLoss.HPO3);
        }

        if(ModificationDB.getInstance().getModificationString(objPeptide).toLowerCase().contains("oxidation of m")){
            specificAnnotationPreferences.addNeutralLoss(NeutralLoss.CH4OS);
            //annotationSettings.addNeutralLoss(NeutralLoss.CH4OS);
        }

        annotationSettings.setIntensityThresholdType(AnnotationParameters.IntensityThresholdType.percentile);

        ModificationParameters modificationParameters = new ModificationParameters();
        SequenceMatchingParameters sequenceMatchingParameters = new SequenceMatchingParameters();
        JSequenceProvider jSequenceProvider = new JSequenceProvider();
        IonMatch []matches = peptideSpectrumAnnotator.getSpectrumAnnotation(annotationSettings,
                specificAnnotationPreferences,
                "",
                "",
                spectrum,
                objPeptide,
                modificationParameters,
                jSequenceProvider,
                sequenceMatchingParameters);


        if ((matches == null || matches.length <=0 ) && PeptideSearchMT.debug) {
            System.err.println("No ions matched!");
            return (null);
        }

        StringBuilder mz_all = new StringBuilder();
        StringBuilder int_all = new StringBuilder();
        StringBuilder mz_match = new StringBuilder();
        StringBuilder int_match = new StringBuilder();
        StringBuilder m_label = new StringBuilder();

        HashMap<Double,Double> validPeakMz = new HashMap<>();

        double[] mzArray = spectrum.mz;
        double[] intensityArray = spectrum.intensity;
        for(int i=0;i<mzArray.length;i++){
            validPeakMz.put(mzArray[i],intensityArray[i]);
            if(i==0){
                mz_all.append(mzArray[i]);
                int_all.append(intensityArray[i]);
            }else{
                mz_all.append(";"+mzArray[i]);
                int_all.append(";"+intensityArray[i]);
            }
        }

        ArrayList<IonMatch> ionMatches = new ArrayList<>();
        for(IonMatch match: matches){
            ionMatches.add(match);
        }

        Iterator<IonMatch> pListIterator = ionMatches.iterator();
        HashSet<Double> ionHitSet = new HashSet<>();
        while (pListIterator.hasNext()) {
            IonMatch ionMatch = pListIterator.next();
            try {
                if (!validPeakMz.containsKey(ionMatch.peakMz) || ionHitSet.contains(ionMatch.peakMz) || (ionMatch.ion.hasNeutralLosses() && ionMatch.ion.getNeutralLosses().length >= 2)) {
                    pListIterator.remove();

                } else {
                    ionHitSet.add(ionMatch.peakMz);
                }
            }catch (NullPointerException nullPointerException){
                nullPointerException.printStackTrace();
                System.exit(1);
            }
        }



        for(int i=0;i<ionMatches.size();i++){
            IonMatch ionMatch = ionMatches.get(i);
            String label = ionMatch.getPeakAnnotation();

            if(i==0){
                mz_match.append(ionMatch.peakMz);
                int_match.append(ionMatch.peakIntensity);
                m_label.append(label);
            }else{
                mz_match.append(";"+ionMatch.peakMz);
                int_match.append(";"+ionMatch.peakIntensity);
                m_label.append(";"+label);
            }

        }

        //IonFactory.getInstance().getDefaultNeutralLosses().clear();

        StringBuilder outline = new StringBuilder();
        outline.append(objPeptide.getSequence()).append("\t");
        outline.append(ModificationDB.getInstance().getModificationString(objPeptide)).append("\t");
        outline.append(spectrum.getSpectrumTitle()).append("\t");
        outline.append(JPeptide.getMass(objPeptide)).append("\t");
        outline.append(spectrum.getPrecursor().mz).append("\t");
        outline.append(spectrum.getPrecursor().possibleCharges[0]).append("\t");
        // previous this is peptide + mod
        outline.append(objPeptide.getSequence()).append("\t");
        outline.append(m_label).append("\t");
        outline.append(mz_match).append("\t");
        outline.append(int_match).append("\t");
        outline.append(mz_all).append("\t");
        outline.append(int_all);

        return(outline.toString());
    }


    public static String getPeakAnnotationHeader(){
        StringBuilder outline = new StringBuilder();
        outline.append("peptide").append("\t");
        outline.append("modification").append("\t");
        outline.append("Query").append("\t");
        outline.append("calc_mr").append("\t");
        outline.append("observed_mz").append("\t");
        outline.append("charge").append("\t");
        outline.append("pepSeq").append("\t");
        outline.append("m_label").append("\t");
        outline.append("m_mz").append("\t");
        outline.append("m_intensity").append("\t");
        outline.append("mz").append("\t");
        outline.append("intensity");

        return(outline.toString());
    }


    public static void plot(String input, String figfile, double limit) throws InterruptedException, IOException {

        InputStream is = PeakAnnotation.class.getResourceAsStream("/main/java/plot/plot_spectrum.r");
        BufferedReader br = new BufferedReader(new InputStreamReader(is));

        String prog_dir = PeakAnnotation.class.getProtectionDomain().getCodeSource().getLocation().getFile();
        File prog_dir_f = new File(prog_dir);
        String sfile = prog_dir_f.getParent();
        sfile = sfile.replaceAll("\\\\", "/");

        String rcode = "arg=c('" + input + "','" + figfile + "','" + limit +"')\n";
        String line = "";
        while ((line = br.readLine()) != null) {
            rcode = rcode + "\n" + line;
        }
        br.close();


        try {

            String rbin= "R";

            Process p = null;
            if(System.getProperty("os.name").equals("Linux")){
                rbin= ""+rbin+"";
                rbin = rbin + " --vanilla --slave";
                p = Runtime.getRuntime().exec(rbin);
            }else if (System.getProperty("os.name").startsWith("Windows")) {
                rbin= "\""+rbin+"\"";
                rbin = rbin + " --vanilla --slave";
                rbin = "cmd /c "+rbin;
                Cloger.getInstance().logger.info(rbin);
                p = Runtime.getRuntime().exec(rbin);
            }else{
                Cloger.getInstance().logger.error("Don't support current System!");
                System.exit(0);
            }

            BufferedReader reader = new BufferedReader(new InputStreamReader(
                    p.getInputStream()));
            BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(
                    p.getOutputStream()));

            BufferedInputStream errin = new BufferedInputStream(p.getErrorStream());
            BufferedReader errbr = new BufferedReader(new InputStreamReader(errin));

            writer.write(rcode + "\n");
            writer.flush();
            writer.close();
            String logline;
            p.waitFor();

            while ((logline = reader.readLine()) != null) {
                Cloger.getInstance().logger.info(logline);
            }


            while ((logline = errbr.readLine()) != null) {
                Cloger.getInstance().logger.info(logline);
            }


        } catch (IOException e) {
            e.printStackTrace();
        }

    }


    public static void plotPSM(String psm_file,String mgf,String outfile) throws IOException {

        ModificationDB.getInstance();


        File mgfFile = new File(mgf);
        WaitingHandler waitingHandler = new WaitingHandlerCLIImpl();
        waitingHandler.setDisplayProgress(false);

        MgfFileIterator mgfFileIterator = new MgfFileIterator(mgfFile, waitingHandler);
        String title;
        Spectrum spectrum = new Spectrum();
        String spectrumTitle;
        HashMap<String,Spectrum> spectrumHashMap = new HashMap<>();
        while ((title = mgfFileIterator.next()) != null) {

            spectrum = mgfFileIterator.getSpectrum();
            spectrum.spectrumTitle = title;
            spectrumHashMap.put(title,spectrum);

        }


        BufferedWriter bWriter = new BufferedWriter(new FileWriter(new File(outfile)));
        bWriter.write(getPeakAnnotationHeader()+"\n");

        BufferedReader bReader = new BufferedReader(new FileReader(new File(psm_file)));

        String headline = bReader.readLine().trim();
        String hd[] = headline.split("\t");
        HashMap<String,Integer> headMap = new HashMap<>();
        for(int i=0;i<hd.length;i++){
            headMap.put(hd[i],i);
        }
        String line;
        while((line = bReader.readLine())!=null){
            String d[] = line.trim().split("\t");
            String peptide = d[headMap.get("peptide")];
            String mod = d[headMap.get("modification")];
            String spectrum_title = d[headMap.get("spectrum_title")];
            spectrum = spectrumHashMap.get(spectrum_title);

            // ArrayList<ModificationMatch> modificationMatches = new ArrayList<>();
            Cloger.getInstance().logger.info(peptide+"\t"+spectrum_title+"\t");
            Peptide ePeptide = new Peptide(peptide);
            if(!mod.equalsIgnoreCase("-")){
                String mods[] = mod.split(";");
                for(String modification : mods){
                    // CLIP_TRAQ_4 of K@13[244.1015]
                    String name = modification.replaceAll("@.*$","");
                    String pos = modification.replaceAll(".*@(\\d+).*$","$1");
                    int modpos = Integer.valueOf(pos);
                    Cloger.getInstance().logger.info(name+"@"+modpos+",");
                    ModificationMatch modificationMatch = new ModificationMatch(name,modpos);
                    // modificationMatches.add(modificationMatch);
                    ePeptide.addVariableModification(modificationMatch);
                }
            }
            Cloger.getInstance().logger.info("");
            //Peptide ePeptide = new Peptide(peptide,modificationMatches);
            String psmMatch_tmp = PeakAnnotation.getPeakAnnotation(ePeptide, spectrum, true);
            bWriter.write(psmMatch_tmp+"\n");

        }

        bReader.close();
        bWriter.close();

    }
}
