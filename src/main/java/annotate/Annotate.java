package main.java.annotate;

import main.java.pg.CParameter;
import main.java.util.Cloger;
import net.sf.jfasta.FASTAElement;
import net.sf.jfasta.FASTAFileReader;
import net.sf.jfasta.impl.FASTAElementIterator;
import net.sf.jfasta.impl.FASTAFileReaderImpl;
import org.apache.commons.cli.*;
import org.apache.commons.lang3.StringUtils;
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;


/**
 * Annotate PepQuery result for taking protein FASTA file as input for searching
 */
public class Annotate {

    public String pepSummaryfile;
    public String vardb;
    public String refdb;
    public String outfile;
    public String pep2profile;
    public String varAnnotationfile;
    public boolean needFilter = false;

    public static void main(String[] args) throws ParseException, IOException {

        Cloger.getInstance().logger.info("Start analysis");
        Cloger.getInstance().logger.info(StringUtils.join(args, " "));

        Options options = new Options();

        //options.addOption("o", true, "Output dir");
        options.addOption("i", true, "PepQuery psm_rank.txt file");
        options.addOption("vardb", true, "Protein sequence file in fasta format");
        options.addOption("refdb", true, "Protein sequence file in fasta format");
        options.addOption("a", true, "Protein annotation file from IPeak");
        options.addOption("m", true, "Peptide to protein mapping file");
        options.addOption("o", true, "Output file");
        options.addOption("f",false,"");

        options.addOption("h", false, "Help");

        CommandLineParser parser = new PosixParser();
        CommandLine cmd = parser.parse(options, args);
        if (cmd.hasOption("h") || cmd.hasOption("help") || args.length == 0) {
            HelpFormatter f = new HelpFormatter();

            // System.out.println(version);
            f.printHelp("Options", options);
            System.exit(0);
        }

        Annotate annotate = new Annotate();
        annotate.pepSummaryfile = cmd.getOptionValue("i");
        annotate.vardb = cmd.getOptionValue("vardb");
        annotate.refdb = cmd.getOptionValue("refdb");
        annotate.pep2profile = cmd.getOptionValue("m");
        annotate.varAnnotationfile = cmd.getOptionValue("a");
        annotate.outfile = cmd.getOptionValue("o");
        if(cmd.hasOption("f")){
            annotate.needFilter = true;
        }else{
            annotate.needFilter = false;
        }
        annotate.run();


    }

    public void run() throws IOException {

        HashMap<String, JVariant> pep2VariantMap = new HashMap<>();

        // read PSM file
        BufferedReader psmReader = new BufferedReader(new FileReader(new File(this.pepSummaryfile)));
        String headLine = psmReader.readLine().trim();
        String h[] = headLine.split("\t");
        HashMap<String,Integer> psmHeadMap = new HashMap<>();
        for(int i=0;i<h.length;i++){
            psmHeadMap.put(h[i],i);
        }

        String line;

        HashMap<String, JVariantPeptide> vpepIdentificationMap = new HashMap<>();

        while((line = psmReader.readLine())!=null){
            line = line.trim();
            String d[] = line.trim().split("\t");
            String pep = d[psmHeadMap.get("peptide")];

            if(this.needFilter){
                int n_db = Integer.parseInt(d[psmHeadMap.get("n_db")]);
                double pvalue = Double.parseDouble(d[psmHeadMap.get("pvalue")]);
                int rank = Integer.parseInt(d[psmHeadMap.get("rank")]);
                int n_ptm = 0;
                if(psmHeadMap.containsKey("n_ptm")){
                    n_ptm = Integer.parseInt(d[psmHeadMap.get("n_ptm")]);
                }
                String peptideSequence = d[psmHeadMap.get("peptide")];
                boolean psm_validation_passed = CParameter.passed_psm_validation(11,pvalue,n_db,rank,n_ptm);
                if(psm_validation_passed){
                    // confident
                }else{
                    continue;
                }
            }

            if(vpepIdentificationMap.containsKey(pep)){
                vpepIdentificationMap.get(pep).npsm = vpepIdentificationMap.get(pep).npsm + 1;
                vpepIdentificationMap.get(pep).pvalue.add(Double.valueOf(d[psmHeadMap.get("pvalue")]));
                vpepIdentificationMap.get(pep).score.add(Double.valueOf(d[psmHeadMap.get("score")]));
            }else{
                JVariantPeptide variantPeptide = new JVariantPeptide();
                variantPeptide.npsm = 1;
                variantPeptide.pvalue.add(Double.valueOf(d[psmHeadMap.get("pvalue")]));
                variantPeptide.score.add(Double.valueOf(d[psmHeadMap.get("score")]));
                vpepIdentificationMap.put(pep,variantPeptide);
            }


            if(pep2VariantMap.containsKey(pep)){

            }else{
                JVariant variant = new JVariant();
                variant.peptideSequence = pep;
                pep2VariantMap.put(pep,variant);
            }
        }


        psmReader.close();

        // Get protein IDs for each peptide
        HashMap<String,String> pep2pro = this.peptideMapping();

        // Get protein sequences
        HashMap<String,String> proSeqMap = new HashMap<>();
        for(String pep: pep2pro.keySet()){
            String acc[] = pep2pro.get(pep).split(";");
            for(String pro: acc){
                proSeqMap.put(pro,"");

            }
        }
        // Get variant protein sequences
        getProSeq(this.vardb,proSeqMap);
        // Get reference protein sequences
        HashMap<String, JVariant> vMap = readVarAnnoInfor();
        HashMap<String,String> refProSeqMap = new HashMap<>();
        for(String varPro: proSeqMap.keySet()){
            refProSeqMap.put(vMap.get(varPro).protein_id,"");
        }
        for(String acc: refProSeqMap.keySet()){
            proSeqMap.put(acc,"");
        }
        getProSeq(this.refdb,proSeqMap);

        // For SNV and INDEL variants
        ArrayList<JVariant> vList = new ArrayList<>();
        for(String pep: pep2VariantMap.keySet()){
            String acc[] = pep2pro.get(pep).split(";");
            for(String varProID: acc){
                JVariant jVariant = new JVariant();
                jVariant.peptideSequence = pep;
                jVariant.vID = varProID;
                jVariant.var_pro_id = varProID;
                jVariant.protein_id = vMap.get(varProID).protein_id;
                jVariant.varLine = vMap.get(varProID).varLine;
                //System.out.println(varProID+"\t"+jVariant.varLine);
                vList.add(jVariant);
            }
        }
        preprocessVariantPeptides(vList,proSeqMap);


        // output
        BufferedWriter pWriter = new BufferedWriter(new FileWriter(new File(this.outfile)));
        pWriter.write("Variant_peptide\tspectra_count\tpvalue\tscore\taa_change_first_pos_in_pep\taa_change_last_pos_in_pep\taa_change_before_pep\taa_change_after_pep\t"+JVariant.varlineColumnNames+"\n");
        for(JVariant vpep: vList){
            if(vpep.valid) {
                StringBuilder oBuilder = new StringBuilder();
                oBuilder.append(vpep.peptideSequence).append("\t");
                oBuilder.append(vpepIdentificationMap.get(vpep.peptideSequence).npsm).append("\t");
                oBuilder.append(StringUtils.join(vpepIdentificationMap.get(vpep.peptideSequence).pvalue,";")).append("\t");
                oBuilder.append(StringUtils.join(vpepIdentificationMap.get(vpep.peptideSequence).score,";")).append("\t");
                oBuilder.append(vpep.aa_change_first_pos_in_pep).append("\t");
                oBuilder.append(vpep.aa_change_last_pos_in_pep).append("\t");
                oBuilder.append(vpep.aa_change_before_pep).append("\t");
                oBuilder.append(vpep.aa_change_after_pep).append("\t");
                oBuilder.append(vpep.varLine).append("\n");
                pWriter.write(oBuilder.toString());
            }
        }
        pWriter.close();


    }

    public HashMap<String, JVariant> readVarAnnoInfor() throws IOException {
        // read annotation file
        BufferedReader pReader = new BufferedReader(new FileReader(new File(this.varAnnotationfile)));
        String headLine = pReader.readLine().trim();
        String h[] = headLine.split("\t");
        HashMap<String,Integer> pHeadMap = new HashMap<>();
        for(int i=0;i<h.length;i++){
            pHeadMap.put(h[i],i);
        }
        JVariant.varlineColumnNames = headLine;
        String line;

        HashMap<String, JVariant> vMap = new HashMap<>();

        while((line = pReader.readLine())!=null){
            line = line.trim();
            String d[] = line.trim().split("\t");
            String proID = d[pHeadMap.get("Protein")];
            String pIDinDB = d[pHeadMap.get("Variant_ID")];
            String geneID = d[pHeadMap.get("Gene")];
            JVariant jVariant = new JVariant();
            jVariant.varLine = line;
            jVariant.vID = pIDinDB;
            jVariant.protein_id = proID;
            jVariant.var_pro_id = pIDinDB;
            vMap.put(jVariant.vID,jVariant);

        }
        pReader.close();
        return(vMap);
    }



    public HashMap<String, String> peptideMapping() throws IOException {
        BufferedReader pReader = new BufferedReader(new FileReader(new File(this.pep2profile)));
        String headLine = pReader.readLine().trim();
        String h[] = headLine.split("\t");
        HashMap<String,Integer> pHeadMap = new HashMap<>();
        for(int i=0;i<h.length;i++){
            pHeadMap.put(h[i],i);
        }
        String line;

        HashMap<String, String> vMap = new HashMap<>();

        while((line = pReader.readLine())!=null){
            line = line.trim();
            String d[] = line.trim().split("\t");
            String pep = d[pHeadMap.get("peptide")];
            String pro = d[pHeadMap.get("protein")];
            vMap.put(pep,pro);

        }
        pReader.close();
        return(vMap);

    }


    public void getProSeq(String dbfile, HashMap<String,String> proSeqMap) throws IOException {
        FASTAFileReader reader = new FASTAFileReaderImpl(new File(dbfile));
        FASTAElementIterator it = reader.getIterator();
        int num = 0;
        while (it.hasNext()) {
            FASTAElement el = it.next();
            el.setLineLength(1);
            String headLine[] = el.getHeader().split("\\s+");
            String proID = headLine[0];
            String proSeq = el.getSequence().toUpperCase();

            proSeq = proSeq.replaceAll("\\*$", "");
            //System.out.println(proID);
            if(proSeqMap.containsKey(proID)){
                proSeqMap.put(proID,proSeq);
                num++;
            }
        }
        Cloger.getInstance().logger.info("Read "+num+" protein sequences from file: "+dbfile);
        reader.close();
    }


    public void preprocessVariantPeptides(ArrayList<JVariant> vList, HashMap<String,String> proSeqMap){
        for(JVariant vpep: vList){
            String proseq_wild = proSeqMap.get(vpep.protein_id);
            String proseq_mut  = proSeqMap.get(vpep.var_pro_id);

            String variant_pep_seq = vpep.peptideSequence;

            int start_mut = proseq_mut.indexOf(variant_pep_seq) + 1;
            int end_mut = start_mut + variant_pep_seq.length() -1;
            vpep.start_pos_pro = start_mut;
            vpep.end_pos_pro = end_mut;

            int end_wild = end_mut;
            if(end_wild > proseq_wild.length()){
                end_wild = proseq_wild.length();
            }

            // the beginning index, inclusive.
            // the ending index, exclusive.
            // Please note the length of variant peptide sequence is possible larger than wild type sequence.
            String wild_pep_seq;
            if(start_mut>proseq_wild.length()){
                wild_pep_seq =  "-";
            }else{
                wild_pep_seq = proseq_wild.substring(start_mut-1,end_wild);
            }
            vpep.knownPeptideSequence = wild_pep_seq;

            String aa_mut[] = variant_pep_seq.split("");
            String aa_wil[] = wild_pep_seq.split("");
            //
            int plen;
            if(wild_pep_seq.length() < variant_pep_seq.length()){
                plen = aa_wil.length;
            }else{
                plen = aa_mut.length;
            }


            int aa_change_first_pos = 0;
            int aa_change_last_pos = 0;
            for(int i=0;i<plen;i++){
                if(!aa_mut[i].equalsIgnoreCase(aa_wil[i])){
                    if(aa_change_first_pos == 0){
                        aa_change_first_pos = i+1;
                    }
                    aa_change_last_pos = i+1;
                }
            }


            if(aa_change_first_pos == 0 && wild_pep_seq.length() < variant_pep_seq.length()){
                aa_change_first_pos = wild_pep_seq.length()+1;
            }

            if(wild_pep_seq.length() < variant_pep_seq.length()){
                aa_change_last_pos = variant_pep_seq.length();
            }

            vpep.aa_change_first_pos_in_pep = aa_change_first_pos;
            vpep.aa_change_last_pos_in_pep  = aa_change_last_pos;

            if(aa_change_first_pos == 0 && aa_change_last_pos == 0){
                Cloger.getInstance().logger.info("Warning: " + vpep.peptideSequence+" digestion problem.");
                vpep.valid = false;
                continue;
            }


            String changed_aas = variant_pep_seq.substring(aa_change_first_pos-1,aa_change_last_pos);
            int last_pos;
            if(aa_change_last_pos > wild_pep_seq.length()){
                last_pos = wild_pep_seq.length();
            }else{
                last_pos = aa_change_last_pos;
            }
            String wild_aas = wild_pep_seq.substring(aa_change_first_pos-1,last_pos);
            vpep.aa_change_before_pep = wild_aas;
            vpep.aa_change_after_pep = changed_aas;

            vpep.aa_change_first_pos_in_pro = vpep.aa_change_first_pos_in_pep + vpep.start_pos_pro - 1;
            vpep.aa_change_last_pos_in_pro  = vpep.aa_change_last_pos_in_pep + vpep.start_pos_pro - 1;

        }

    }
}
