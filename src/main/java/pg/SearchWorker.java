package main.java.pg;


import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import main.java.OpenModificationSearch.ModificationDB;
import main.java.PSMMatch.JPeptide;
import main.java.PSMMatch.ScoreResult;
import main.java.util.Cloger;
import java.util.ArrayList;

import static main.java.pg.PeptideSearchMT.scorePeptide2Spectrum;

/**
 * After scoring PSMs for target peptides, scoring PSMs from reference db and random peptides. (multiple treads)
 */
public final class SearchWorker implements Runnable {


    private PeptideInput peptideInput = null;
    private DatabaseInput databaseInput = null;

    private int index;

    public SearchWorker(PeptideInput pepInput, DatabaseInput dbInput, int ind){
        this.peptideInput = pepInput;
        this.databaseInput = dbInput;
        this.index = ind;
    }

    @Override
    public void run() {

        // peptide form without variable modification
        Cloger.getInstance().logger.debug(Thread.currentThread().getName()+" -> Start "+this.index+": "+peptideInput.peptideSequence+" ...");
        // peptide form with variable modification
        for (JPeptide jPeptide : peptideInput.getPtmIsoforms()) {
            for (int i=0;i < jPeptide.spectraIndexs.size();i++) {
                Spectrum spectrum = SpectraInput.spectraMap.get(jPeptide.spectraIndexs.get(i));
                ScoreResult scoreResult = jPeptide.scores.get(i);
                if(!jPeptide.valid.get(i)){
                    continue;
                }

                String spectrumTitle = spectrum.getSpectrumTitle();

                // output the best match from reference database searching
                ScoreResult bestRefScoreResult = new ScoreResult();
                int k = 0;
                if(databaseInput.ms2peptide.containsKey(spectrumTitle)) {
                    for (JPeptide ePeptide : databaseInput.ms2peptide.get(spectrumTitle)) {
                        //scoreList[k] = scorePeptide2Spectrum(ePeptide,spectrum,itol).score;
                        ScoreResult tmpScoreResult = scorePeptide2Spectrum(ePeptide.peptide, spectrum);
                        tmpScoreResult.jPeptide = ePeptide;
                        tmpScoreResult.spectrum_title = spectrumTitle;
                        double tmpScore = tmpScoreResult.score;

                        //double delta_ppm = (JPeptide.getMass(ePeptide) - spectrum.getPrecursor().getMass(spectrum.getPrecursor().possibleCharges[0]))/JPeptide.getMass(ePeptide)*1.0e6;
                        //System.out.println(spectrumTitle+"\t"+ePeptide.getSequence()+"\t"+ModificationDB.getInstance().getModificationString(ePeptide)+"\t"+tmpScore+"\t"+delta_ppm);

                        // save the best matching from the reference database searching for the target spectrum
                        if(bestRefScoreResult.score < tmpScore){
                            bestRefScoreResult = tmpScoreResult;
                        }

                        //scoreList.add(tmpScore);
                        if (scoreResult.score <= tmpScore) {
                            scoreResult.nBetter++;

                            if (tmpScore >= 20) {
                                String dout = get_psm_outline(tmpScoreResult);
                                this.peptideInput.addOutAdditionalLines(dout);
                            }
                            if(CParameter.fast_model){
                                break;
                            }
                        }

                        if(CParameter.scoreMethod == 0){
                            double tmpScore2 = tmpScoreResult.score2;
                            if (scoreResult.score2 <= tmpScore2) {
                                scoreResult.nBetter2++;
                            }
                        }

                        k++;
                    }
                }

                // save the best matching from the reference database searching for the target spectrum
                if(bestRefScoreResult.score > 0){
                    String dout = get_psm_outline(bestRefScoreResult);
                    if(!this.peptideInput.getOutAdditionalLines().contains(dout)){
                        this.peptideInput.addOutAdditionalLines(dout);
                    }

                }


                double pvalue = 100;
                double pvalue2 = 100;
                int nRand = -1;
                int nRand2 = -1;

                ArrayList<Double> randomScoreList1 = new ArrayList<>(10000);
                ArrayList<Double> randomScoreList2 = new ArrayList<>(10000);

                if(scoreResult.nBetter==0 && scoreResult.nBetter2==0){
                    nRand = 0;
                    nRand2 = 0;
                    if(peptideInput.randomPeptides.size()<=0){
                        peptideInput.setRandomPeptides(PeptideInput.generateRandomPeptidesFast(peptideInput.peptideSequence,CParameter.nRandomPeptides,true));
                    }
                    for (String randomPeptideSeq : peptideInput.randomPeptides) {
                        Peptide randomPeptide = PeptideInput.getModificationPeptide(jPeptide.peptide, randomPeptideSeq, CParameter.fixMods);
                        ScoreResult tmpScoreResult = scorePeptide2Spectrum(randomPeptide, spectrum);
                        double randomScore = tmpScoreResult.score;

                        randomScoreList1.add(randomScore);

                        if (randomScore >= scoreResult.score) {
                            nRand++;
                        }

                        if(CParameter.scoreMethod == 0){
                            double randomScore2 = tmpScoreResult.score2;
                            randomScoreList2.add(randomScore2);
                            if (randomScore2 >= scoreResult.score2) {
                                nRand2++;
                            }

                        }

                    }
                    pvalue = 1.0 * (nRand + 1) / (peptideInput.randomPeptides.size() + 1);

                    if(CParameter.scoreMethod ==0){
                        pvalue2 = 1.0 * (nRand2 + 1) / (peptideInput.randomPeptides.size() + 1);
                    }


                }



                int sumRand = -1;
                if(nRand != -1){
                    sumRand = peptideInput.randomPeptides.size();
                }


                StringBuilder outBuilder = new StringBuilder();
                //double tolPpm = (1.0e6) * (jPeptide.getMass() - spectrum.getPrecursor().getMass(spectrum.getPrecursor().possibleCharges[0])) / jPeptide.getMass();
                ArrayList<Double> tol_res = CParameter.get_mass_error(spectrum.getPrecursor().getMass(spectrum.getPrecursor().possibleCharges[0]), jPeptide.getMass());

                /**
                if(tol_res.size()<=0){
                    Cloger.getInstance().logger.error(jPeptide.peptide.getSequence());
                    Cloger.getInstance().logger.error(spectrum.getSpectrumTitle());
                    Cloger.getInstance().logger.error(spectrum.getPrecursor().getMass(spectrum.getPrecursor().possibleCharges[0]));
                    Cloger.getInstance().logger.error(jPeptide.getMass());
                }**/

                if(CParameter.scoreMethod==0){
                    pvalue = Math.max(pvalue,pvalue2);
                }

                int nBetter = Math.max(scoreResult.nBetter,scoreResult.nBetter2);
                outBuilder.append(jPeptide.peptide.getSequence())
                        .append("\t")
                        //.append(jPeptide.peptide.getSequenceWithLowerCasePtms())
                        .append(ModificationDB.getInstance().getModificationString(jPeptide.peptide))
                        .append("\t")
                        .append(jPeptide.spectraIndexs.size())
                        .append("\t")
                        .append(spectrum.getSpectrumTitle())
                        .append("\t")
                        .append(spectrum.getPrecursor().possibleCharges[0])
                        .append("\t")
                        .append(spectrum.getPrecursor().getMass(spectrum.getPrecursor().possibleCharges[0]))
                        .append("\t")
                        .append(tol_res.get(0))
                        .append("\t")
                        .append(tol_res.get(1))
                        .append("\t")
                        .append(tol_res.get(2))
                        .append("\t")
                        .append(jPeptide.getMass())
                        .append("\t")
                        .append(spectrum.getPrecursor().mz)
                        .append("\t")
                        .append(scoreResult.score)
                        .append("\t")
                        .append(nBetter)
                        .append("\t")
                        .append(k)
                        .append("\t")
                        .append(nRand)
                        .append("\t")
                        .append(sumRand)
                        .append("\t")
                        .append(pvalue);


                //System.out.println(outBuilder.toString());
                this.peptideInput.addOutputLine(outBuilder.toString());


            }

            jPeptide.scores = null;

        }

        this.peptideInput.randomPeptides = null;

        Cloger.getInstance().logger.debug(Thread.currentThread().getName()+" -> Finished "+this.index+": "+peptideInput.peptideSequence+" .");




    }

    /**
     * Export information for a given PSM. This is used to export PSMs from reference database matching.
     * @param scoreResult A ScoreResult object
     * @return A string contains PSM information.
     */
    public static String get_psm_outline(ScoreResult scoreResult){
        Peptide ePeptide = scoreResult.jPeptide.peptide;
        String spectrum_title = scoreResult.spectrum_title;
        Spectrum spectrum = SpectraInput.getSpectrum(spectrum_title);
        //double pep_mass = JPeptide.getMass(ePeptide);
        double pep_mass = scoreResult.jPeptide.getMass();
        // double pep_mass = ePeptide.getMass();
        double pre_mass = spectrum.getPrecursor().getMass(spectrum.precursor.possibleCharges[0]);
        ArrayList<Double> tol_res = CParameter.get_mass_error(pre_mass,pep_mass);

        if(tol_res.size()<=0){
            Cloger.getInstance().logger.error("Mass tol error: "+spectrum.getSpectrumTitle() + "\t" +
                    ePeptide.getSequence() + "\t" +
                    //ePeptide.getSequenceWithLowerCasePtms() + "\t" +
                    ModificationDB.getInstance().getModificationString(ePeptide)+ "\t" +
                    pre_mass + "\t" +
                    pep_mass + "\n");
            System.exit(1);
        }

        String dout = spectrum.getSpectrumTitle() + "\t" +
                ePeptide.getSequence() + "\t" +
                //ePeptide.getSequenceWithLowerCasePtms() + "\t" +
                ModificationDB.getInstance().getModificationString(ePeptide)+ "\t" +
                pre_mass + "\t" +
                pep_mass + "\t" +
                tol_res.get(0) + "\t" +
                tol_res.get(1) + "\t" +
                tol_res.get(2) + "\t" +
                scoreResult.score;
        return dout;
    }
}
