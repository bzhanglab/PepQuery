package main.java.pg;

import main.java.PSMMatch.JPeptide;

import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

public class StatisticScoreCalc {

    public ArrayList<PeptideInput> peptideInputs = new ArrayList<>();

    public StatisticScoreCalc(ArrayList<PeptideInput> peptideInputs){
        this.peptideInputs = peptideInputs;
    }


    public void run() throws InterruptedException {
        ExecutorService fixedThreadPool = Executors.newFixedThreadPool(CParameter.cpu);

        for(PeptideInput peptideInput: peptideInputs) {
            // PeptideInput peptideInput = peptideInputs.get(0);
            boolean require_p_value = false;
            for (JPeptide jPeptide : peptideInput.getPtmIsoforms()) {
                for (int i = 0; i < jPeptide.spectraIndexs.size(); i++) {
                    if (jPeptide.valid.get(i)) {
                        require_p_value = true;
                        break;
                    }
                }
                if(require_p_value){
                    break;
                }
            }
            if(require_p_value) {
                peptideInput.setRandomPeptides(PeptideInput.generateRandomPeptidesFast(peptideInput.getPtmIsoforms().get(0).peptide.getSequence(), CParameter.nRandomPeptides, true));
            }
        }

        fixedThreadPool.shutdown();
        fixedThreadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);

        fixedThreadPool = Executors.newFixedThreadPool(CParameter.cpu);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for(PeptideInput peptideInput: peptideInputs) {
            // PeptideInput peptideInput = peptideInputs.get(0);
            //peptideInput.setRandomPeptides(PeptideInput.generateRandomPeptidesFast(peptideInput.getPtmIsoforms().get(0).peptide.getSequence(), CParameter.nRandomPeptides, true));
            for (JPeptide jPeptide : peptideInput.getPtmIsoforms()) {
                //for (int i=0;i<jPeptide.mSnSpectrums.size();i++) {
                for (int i = 0; i < jPeptide.spectraIndexs.size(); i++) {
                    // public PsmScoringWorker (PeptideInput pi, DatabaseInput di, int index, double itol, int method,String fixModOptionValue, JPeptide jPeptide){
                    if (jPeptide.valid.get(i)){
                        //fixedThreadPool.execute(new PsmScoringWorker(peptideInput, databaseInput, i, jPeptide));
                        fixedThreadPool.execute(new StatisticScoreWorker(peptideInput, jPeptide, i));
                    }
                }
            }
        }

        fixedThreadPool.shutdown();
        fixedThreadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);
    }
}
