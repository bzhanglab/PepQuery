package main.java.pg;

import com.compomics.util.experiment.biology.modifications.Modification;
import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.identification.protein_sequences.SingleProteinSequenceProvider;
import com.compomics.util.experiment.identification.utils.ModificationUtils;
import com.compomics.util.experiment.io.biology.protein.SequenceProvider;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import main.java.OpenModificationSearch.ModificationDB;
import main.java.util.MSDataSet;
import org.apache.commons.io.IOUtils;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;


public final class ModificationGear {

    private static ModificationGear instance = null;
    private static SequenceProvider sequenceProvider = new SingleProteinSequenceProvider();
    private static SequenceMatchingParameters modificationsSequenceMatchingParameters = SequenceMatchingParameters.getDefaultSequenceMatching();


    private HashMap<Integer,String> id2ptmname = new HashMap<Integer, String>();


    private ModificationGear(){
        ModificationFactory ptmFactory = ModificationFactory.getInstance();
        //
        //ArrayList<String> ptmNames = ptmFactory.getModifications();
        ArrayList<String> ptmNames = ptmFactory.getDefaultModificationsOrdered();
        ArrayList<String> top_mods = load_top_modifications();
        ptmNames.removeAll(top_mods);
        top_mods.addAll(ptmNames);
        int i=0;
        HashSet<String> ptm_name_list = new HashSet<>();
        for(String name : top_mods){
            Modification ptm = ptmFactory.getModification(name);
            if(!ptm.getCategory().name().contains("Nucleotide_Substitution")) {
                i++;
                id2ptmname.put(i, name);
                ptm_name_list.add(name);
            }
        }


        // the alphabetically ordered names of the user defined modifications
        ModificationDB.getInstance();
        ArrayList<String> usrPTMs = ptmFactory.getUserModificationsOrdered();
        for (String name: usrPTMs){
            if(ptm_name_list.contains(name)){
                continue;
            }
            Modification ptm = ptmFactory.getModification(name);
            if(!ptm.getCategory().name().contains("Nucleotide_Substitution")) {
                i++;
                id2ptmname.put(i, name);
            }
        }


    }

    public static ModificationGear getInstance() {
        if (instance == null) {
            instance = new ModificationGear();
        }
        return instance;
    }

    /**
     *
     * @return A list of PTM names
     */
    public ArrayList<String> load_top_modifications(){

        InputStream inputStream = MSDataSet.class.getResourceAsStream("/main/resources/top_modifications.tsv");
        ArrayList<String> ptm_names = new ArrayList<>();
        try {
            List<String> mods = IOUtils.readLines(inputStream, StandardCharsets.UTF_8);
            String[] head = mods.get(0).split("\t");
            HashMap<String,Integer> column2index = new HashMap<>();
            for(int i=0;i<head.length;i++){
                column2index.put(head[i],i);
            }
            for(int i=1;i<mods.size();i++){
                String[] line = mods.get(i).split("\t");
                ptm_names.add(line[column2index.get("mod_name")]);
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        return ptm_names;
    }

    public Modification getPTM(int index){
        ModificationFactory ptmFactory = ModificationFactory.getInstance();
        return(ptmFactory.getModification(this.id2ptmname.get(index)));
    }


    public String getPTMname(int index){
        ModificationFactory ptmFactory = ModificationFactory.getInstance();
        return(this.id2ptmname.get(index));
    }

    /**
     * Get the ptms according to ptm index (1,2,3)
     * @param mod modification index, such as 1,2,3
     * @return
     */
    public ArrayList<Modification> getPTMs(String mod){
        ArrayList<Modification> ptms = new ArrayList<Modification>();
        ModificationFactory ptmFactory = ModificationFactory.getInstance();
        String[] d = mod.split(",");
        for(int i=0;i<d.length;i++){
            ptms.add(ptmFactory.getModification(this.id2ptmname.get(Integer.valueOf(d[i]))));
        }
        return(ptms);
    }

    /**
     * This function could be used to get fixed modification for a peptide
     * @param pepObj
     * @param modifications
     * @return
     */
    public HashMap<Integer, Modification> getPossibleModificationSites(Peptide pepObj, ArrayList<Modification> modifications){
        HashMap<Integer,Modification> position2ptm = new HashMap<>();
        for (Modification ptm : modifications) {
            //ArrayList<Integer> possibleSites = pepObj.getPotentialModificationSitesNoCombination(ptm, pepObj.getSequence(), 0);
            int pSites[] = ModificationUtils.getPossibleModificationSites(pepObj, ptm, sequenceProvider, modificationsSequenceMatchingParameters);

            if(pSites.length >= 1) {
                ArrayList<Integer> possibleSites = new ArrayList<>();
                for (int i : pSites) {
                    possibleSites.add(i);
                }
                for (Integer i : possibleSites) {
                    position2ptm.put(i, ptm);
                }
            }
        }
        return position2ptm;
    }

    public void printPTM(){
        ModificationFactory ptmFactory = ModificationFactory.getInstance();
        System.out.println("mod_id\tmod_name\tmod_mass\tmod_type\tmod_category\tunimod_accession");
        for(int ptm_id : id2ptmname.keySet()){
            Modification ptm = ptmFactory.getModification(id2ptmname.get(ptm_id)) ;
            String unimod_acc = "-";
            if(ptm.getUnimodCvTerm()!=null){
                unimod_acc = ptm.getUnimodCvTerm().getAccession();
            }
            System.out.println(ptm_id + "\t" + ptm.getName() + "\t" + ptm.getMass() + "\t" + ptm.getModificationType()+"\t"+ptm.getCategory().name()+"\t"+unimod_acc);

        }
    }
}
