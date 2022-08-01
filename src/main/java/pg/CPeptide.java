package main.java.pg;

import java.util.ArrayList;

/**
 * Peptide class for peptides from DNA/protein/VCF information
 */
public class CPeptide {

    public String peptideSequence;
    public int position;
    public int missedClevageSites = -1;
    public int nhit = 0;

    /**
     * When this peptide is from a protein that is translated from a DNA sequence,
     * then the frame will be recorded in this variable. "-1" means this peptide is
     * not from an DNA sequence.
     */
    public int frame = -1;

    /**
     * protein ID
     */
    public ArrayList<String> proteinID = new ArrayList<>();

}
