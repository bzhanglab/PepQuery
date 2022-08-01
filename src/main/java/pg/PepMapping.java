package main.java.pg;


import com.compomics.util.experiment.identification.protein_inference.PeptideProteinMapping;
import com.compomics.util.experiment.identification.protein_inference.fm_index.FMIndex;
import com.compomics.util.experiment.identification.protein_inference.FastaMapper;

import com.compomics.util.gui.waiting.waitinghandlers.WaitingHandlerCLIImpl;

import java.io.*;
import java.util.ArrayList;

import com.compomics.util.parameters.identification.advanced.PeptideVariantsParameters;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import com.compomics.util.parameters.identification.search.ModificationParameters;
import com.compomics.util.parameters.identification.search.SearchParameters;
import main.java.util.Cloger;

public class PepMapping {

    public static void main(String[] args){
        String db = args[0];
        String pep = args[1];

        PepMapping pepMapping = PepMapping.getInstance();
        pepMapping.loadDB(db);
        //pepMapping.close();



    }



    private static PepMapping instance = null;
    private String database;

    public static PepMapping getInstance() {
        if (instance == null) {
            instance = new PepMapping();
        }
        return instance;
    }

    public void loadDB(String db){
        this.database = db;
        this.init();
    }

    /**
    public void close(){
        ((FMIndex)peptideMapper).close();
    }
     **/


    private FastaMapper peptideMapper = null;
    private SequenceMatchingParameters sequenceMatchingParameters = new SequenceMatchingParameters();

    private void init(){

        Cloger.getInstance().logger.info("Load db file: "+database);
        File fastaFile = new File(this.database);

        WaitingHandlerCLIImpl waitingHandlerCLIImpl = new WaitingHandlerCLIImpl();

        SearchParameters searchParameters = null;
        PeptideVariantsParameters peptideVariantsPreferences = PeptideVariantsParameters.getNoVariantPreferences();

        searchParameters = new SearchParameters();
        searchParameters.setModificationParameters(new ModificationParameters());
        searchParameters.setFragmentIonAccuracy(0.02);
        searchParameters.setFragmentAccuracyType(SearchParameters.MassAccuracyType.DA);
        sequenceMatchingParameters = new SequenceMatchingParameters();
        sequenceMatchingParameters.setSequenceMatchingType(SequenceMatchingParameters.MatchingType.indistiguishableAminoAcids);
        sequenceMatchingParameters.setLimitX(0.25);

        Cloger.getInstance().logger.info("Start indexing fasta file");
        long startTimeIndex = System.nanoTime();
        try {
            peptideMapper = new FMIndex(fastaFile, null, waitingHandlerCLIImpl, false, peptideVariantsPreferences, searchParameters);
        } catch (IOException e) {
            Cloger.getInstance().logger.error("Error: cound not index the fasta file");
            e.printStackTrace();
            System.exit(-1);
        }
        double diffTimeIndex = System.nanoTime() - startTimeIndex;
        Cloger.getInstance().logger.info("Indexing took " + (diffTimeIndex / 1e9) + " seconds and consumes " + (((float) ((FMIndex) peptideMapper).getAllocatedBytes()) / 1e6) + " MB");

    }

    /**
     *
     * @param peptideSequence peptide sequence
     * @return
     */
    public boolean hasProtein(String peptideSequence){

        boolean match = false;

        ArrayList<PeptideProteinMapping> peptideProteinMappings = peptideMapper.getProteinMapping(peptideSequence, sequenceMatchingParameters);
        if(peptideProteinMappings.size() >0 ){
            match = true;
        }

        return match;
    }


    /**
     * Retrieve protein information
     * @param peptideSequence peptide sequence
     * @return
     */
    public ArrayList<PeptideProteinMapping> retrieveProtein(String peptideSequence){

        ArrayList<PeptideProteinMapping> peptideProteinMappings = peptideMapper.getProteinMapping(peptideSequence, sequenceMatchingParameters);
        return peptideProteinMappings;
    }
}
