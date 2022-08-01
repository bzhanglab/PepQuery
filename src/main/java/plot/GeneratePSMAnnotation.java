package main.java.plot;

import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import main.java.PSMMatch.JPeptide;
import main.java.pg.CParameter;
import main.java.pg.PeptideInput;
import main.java.pg.RankPSM;
import main.java.pg.SpectraInput;
import main.java.util.Cloger;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GeneratePSMAnnotation {

    public void run(ArrayList<PeptideInput> peptideInputs) throws InterruptedException, IOException {

        ExecutorService fixedThreadPoolScore = Executors.newFixedThreadPool(CParameter.cpu);

        ConcurrentHashMap<Integer,String> annoMap = new ConcurrentHashMap<>(100);
        int i=0;
        for(PeptideInput peptideInput : peptideInputs) {
            for (JPeptide jPeptide : peptideInput.getPtmIsoforms()) {
                for (String title : jPeptide.spectraIndexs) {
                    Spectrum spectrum = SpectraInput.spectraMap.get(title);
                    Peptide peptide = jPeptide.peptide;
                    i = i + 1;
                    annoMap.put(i,"-");
                    fixedThreadPoolScore.execute(new PSMAnnoWorker(peptide, spectrum, true, annoMap, i));

                }
            }

        }

        fixedThreadPoolScore.shutdown();
        fixedThreadPoolScore.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);

        String psm_anno_file = CParameter.outdir + "/psm_annotation.txt";
        File P = new File(psm_anno_file);

        BufferedWriter annoWriter = new BufferedWriter(new FileWriter(psm_anno_file,true));
        if(!P.isFile()){
            annoWriter.write(PeakAnnotation.getPeakAnnotationHeader()+"\n");
        }
        for(int k:annoMap.keySet()){
            if (!annoMap.get(k).equals("-")) {
                annoWriter.write(annoMap.get(k) + "\n");
            }
        }
        annoWriter.close();
        annoMap.clear();
    }

    public void annotateBestPsmFromModSearch(String ptm_file) throws IOException, InterruptedException {

        Cloger.getInstance().logger.info("Annotate best PSMs from unrestricted modification searching ...");

        File ptm_F = new File(ptm_file);
        if(!ptm_F.exists()){
            Cloger.getInstance().logger.info("File:"+ptm_file+" doesn't exist!");
            return;
        }

        // ptm.txt
        // spectrum_title  peptide charge  exp_mass        pep_mass        modification    score
        // read PSM data
        BufferedReader psmReader = new BufferedReader(new FileReader(ptm_file));
        String headline = psmReader.readLine().trim();
        String[] h = headline.split("\t");
        HashMap<String,Integer> hMap = new HashMap<>();
        for(int i=0;i<h.length;i++){
            hMap.put(h[i],i);
        }
        String line;
        HashMap<String, RankPSM> bestPsmMap = new HashMap<>();

        while((line = psmReader.readLine())!=null){
            line = line.trim();
            String[] d = line.split("\t");
            String spectrum_title = d[hMap.get("spectrum_title")];
            double score = Double.parseDouble(d[hMap.get("score")]);

            if(bestPsmMap.containsKey(spectrum_title)){
                if(bestPsmMap.get(spectrum_title).score < score){
                    bestPsmMap.get(spectrum_title).score = score;
                    bestPsmMap.get(spectrum_title).line = line;
                }
            }else{
                RankPSM rankPSM = new RankPSM();
                rankPSM.score = score;
                rankPSM.line = line;
                bestPsmMap.put(spectrum_title,rankPSM);
            }


        }

        psmReader.close();

        ExecutorService fixedThreadPoolScore = Executors.newFixedThreadPool(CParameter.cpu);

        ConcurrentHashMap<Integer,String> annoMap = new ConcurrentHashMap<>(100);
        int i=0;

        Pattern pattern = Pattern.compile("^(.*)@(\\d+)");
        for(String spectrum_title: bestPsmMap.keySet()){
            String[] d = bestPsmMap.get(spectrum_title).line.split("\t");
            String peptide_seq = d[hMap.get("peptide")];
            String modification = d[hMap.get("modification")];
            //ArrayList<ModificationMatch> modificationMatches = new ArrayList<>();
            Peptide peptide = new Peptide(peptide_seq);
            if(!modification.startsWith("-")){
                String[] mods = modification.split(";");
                // TMT2plex of S@4[225.1558]
                for(String mod : mods){
                    Matcher matcher = pattern.matcher(mod);
                    while (matcher.find()) {
                        String ptm_name = matcher.group(1);
                        int pos = Integer.parseInt(matcher.group(2));
                        // PTM ptm = ptmFactory.getPTM(ptm_name);
                        ModificationMatch mm = new ModificationMatch(ptm_name, pos);
                        //modificationMatches.add(mm);
                        peptide.addVariableModification(mm);
                    }
                }
            }

            Spectrum spectrum = SpectraInput.spectraMap.get(spectrum_title);
            //String psmMatch_tmp = PeakAnnotation.getPeakAnnotation(peptide, spectrum, true);
            i++;
            annoMap.put(i,"-");
            fixedThreadPoolScore.execute(new PSMAnnoWorker(peptide, spectrum, true, annoMap, i));

        }

        fixedThreadPoolScore.shutdown();
        fixedThreadPoolScore.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);

        String psm_anno_file = CParameter.outdir + "/psm_annotation.txt";
        File P = new File(psm_anno_file);

        BufferedWriter annoWriter = new BufferedWriter(new FileWriter(psm_anno_file,true));
        if(!P.isFile()){
            annoWriter.write(PeakAnnotation.getPeakAnnotationHeader()+"\n");
        }
        for(int k:annoMap.keySet()){
            if (!annoMap.get(k).equals("-")) {
                annoWriter.write(annoMap.get(k) + "\n");
            }
        }
        annoWriter.close();
        annoMap.clear();

        Cloger.getInstance().logger.info("Generate annotation data for best matches from unrestricted modification searching:"+bestPsmMap.size());

    }

    public class PSMAnnoWorker implements  Runnable {

        private Peptide peptide;
        private Spectrum spectrum;
        boolean lossWaterNH3;

        ConcurrentHashMap<Integer,String> annoMap;

        int i;

        public PSMAnnoWorker(Peptide peptide, Spectrum spectrum, boolean lossWaterNH3, ConcurrentHashMap<Integer,String> annoMap, int i){
            this.peptide = peptide;
            this.spectrum = spectrum;
            this.lossWaterNH3 = lossWaterNH3;
            this.annoMap = annoMap;
            this.i = i;
        }

        @Override
        public void run() {
            String psmMatch_tmp = PeakAnnotation.getPeakAnnotation(peptide, spectrum, true);
            if(psmMatch_tmp!=null){
                this.annoMap.put(i,psmMatch_tmp);
            }
        }
    }


}
