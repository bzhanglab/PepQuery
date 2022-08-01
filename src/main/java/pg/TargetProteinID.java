package main.java.pg;


import main.java.util.Cloger;
import net.sf.jfasta.FASTAElement;
import net.sf.jfasta.FASTAFileReader;
import net.sf.jfasta.impl.FASTAElementIterator;
import net.sf.jfasta.impl.FASTAFileReaderImpl;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * This class is used to handle target protein identification
 */
public class TargetProteinID {


    /**
     * The default value is "-". This is used for checking if there is sequence assigned to it.
     * If not, the value should be "-".
     */
    private String targetProteinSequence = "-";
    private String targetProteinAcc;
    private String refdb;

    /**
     * In target protein identification mode, only identity the decoy version of selected protein.
     * This is used to evaluate the FDR for large scale target protein identification.
     */
    public boolean only_identity_decoy_version = false;


    /**
     *
     * @param db database
     * @param accession A single protein ID or multiple protein IDs separated by ","
     * @param outdir output dir
     * @throws IOException
     */
    public void prepareDB(String db, String accession, String outdir) throws IOException {
        //

        this.refdb = outdir + "/refdb.fasta";
        BufferedWriter dbWriter = new BufferedWriter(new FileWriter(new File(this.refdb)));
        this.targetProteinAcc = accession;

        HashMap<String,Integer> targeted_protein_IDs = new HashMap<>();
        String [] protein_acc_list = accession.split(",");
        for(String acc : protein_acc_list){
            targeted_protein_IDs.put(acc, 0);
        }

        File dbFile = new File(db);
        FASTAFileReader reader = new FASTAFileReaderImpl(dbFile);
        FASTAElementIterator it = reader.getIterator();
        int num = 0;
        while (it.hasNext()) {
            FASTAElement el = it.next();
            el.setLineLength(1);
            String headLine[] = el.getHeader().split("\\s+");
            String proID = headLine[0];
            String proSeq = el.getSequence().toUpperCase();
            //if(proID.equals(accession)){
            if(targeted_protein_IDs.containsKey(proID)){
                targeted_protein_IDs.put(proID,1);
                if(this.only_identity_decoy_version){
                    if(targeted_protein_IDs.size()>=2){
                        Cloger.getInstance().logger.error("Decoy protein searching doesn't support multiple proteins as input!");
                        System.exit(1);
                    }
                    // identity decoy protein
                    StringBuilder stringBuilder = new StringBuilder();
                    stringBuilder.append(proSeq);
                    if(CParameter.decoy_method.equalsIgnoreCase("reverse")) {
                        Cloger.getInstance().logger.info("Generate decoy protein by reverse method");
                        this.targetProteinSequence = stringBuilder.reverse().toString();
                    }else if(CParameter.decoy_method.equalsIgnoreCase("random")){
                        Cloger.getInstance().logger.info("Generate decoy protein by random method");
                        List<String> aalist = Arrays.asList(proSeq.split(""));
                        Collections.shuffle(aalist);
                        this.targetProteinSequence = String.join("",aalist);
                    }else{
                        Cloger.getInstance().logger.error("Generate decoy protein by a method not supported yet:"+CParameter.decoy_method);
                        System.exit(1);
                    }


                    dbWriter.write(">"+el.getHeader()+"\n");
                    dbWriter.write(proSeq+"\n");

                }else {
                    // identity target protein
                    if(this.targetProteinSequence.equalsIgnoreCase("-")){
                        this.targetProteinSequence = proSeq;
                    }else{
                        // multiple protein IDs as input
                        this.targetProteinSequence = this.targetProteinSequence+","+proSeq;
                    }

                    // Cloger.getInstance().logger.info("Try to identify protein: ("+accession+") "+this.targetProteinSequence);
                }
            }else{
                dbWriter.write(">"+el.getHeader().split("\\s")[0]+"\n");
                dbWriter.write(proSeq+"\n");
            }
            num++;
            // System.out.println(proID);
        }
        reader.close();
        dbWriter.close();

        int n_targeted_protein_valid = 0;
        int n_targeted_protein_invalid = 0;
        for(String pro: targeted_protein_IDs.keySet()){
            if(targeted_protein_IDs.get(pro) < 1){
                n_targeted_protein_invalid++;
                Cloger.getInstance().logger.warn("Target protein is not present in database ("+db+"): "+pro+", ignored!");
                //System.exit(1);
            }else{
                n_targeted_protein_valid++;
            }
        }
        Cloger.getInstance().logger.info("Input protein(s):"+targeted_protein_IDs.size()+", valid protein(s):"+n_targeted_protein_valid+", invalid protein(s):"+n_targeted_protein_invalid);


        if(!this.targetProteinSequence.equalsIgnoreCase("-")){
            Cloger.getInstance().logger.info("Try to identity protein: (" + accession + ") " + this.targetProteinSequence);
        }else{
            Cloger.getInstance().logger.error("No valid protein found in database:" + db);
            System.exit(1);
        }


        Cloger.getInstance().logger.info("Proteins: "+num);
    }

    public String getRefdb() {
        return refdb;
    }

    public String getTargetProteinSequence() {
        return targetProteinSequence;
    }

}
