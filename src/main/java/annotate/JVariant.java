package main.java.annotate;


public class JVariant {


    public String vID;
    public String line;

    public static String headLine;

    /**
     * start position in protein, it's 1-based.
     */
    public int start_pos_pro;


    /**
     * end position in protein, it's 1-based.
     */
    public int end_pos_pro;


    /**
     * This is the wild type protein ID
     */
    public String protein_id;

    /**
     * This is variant protein ID which is from customized database
     */
    public String var_pro_id;


    /**
     * Variant type:
     * 1=single amino acid change, from amino acid to amino acid;
     * 2=single amino acid change, from amino acid to stop coden;
     * 3=single amino acid change, from stop coden to amino acid;
     * 4=
     */
    // public int variant_type;


    public String peptideSequence;

    /**
     * The corresponding known peptide sequence in the same position
     */
    public String knownPeptideSequence;

    // public String proteinSequence;


    public int aa_change_first_pos_in_pep;


    public int aa_change_last_pos_in_pep;


    public int aa_change_first_pos_in_pro;


    public int aa_change_last_pos_in_pro;




    public String aa_change_before_pep = "-";
    public String aa_change_after_pep = "-";

    /**
     *
     */
    public String rawLine;

    public String varLine;

    public static String varlineColumnNames;

    public boolean valid = true;

}
