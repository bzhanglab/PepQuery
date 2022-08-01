package main.java.PSMMatch;

import com.compomics.util.experiment.biology.proteins.Peptide;


public class ScoreResult {

    /**
     * The Primary score. This score is used to filter PSMs in searching.
     * In default, this is hyperscore.
     */
    public double score = -1;

    public int nBetter = 0;

    /**
     * The second score. We use two scoring systems.
     * In default, this is MVH score.
     */
    public double score2 = -1;

    /**
     * It's used for score2
     */
    public int nBetter2 = 0;


    // reference database searching related values
    public int n_matched_ref_peptides = 0;
    public int n_matched_ref_peptides_with_better_score = 0;

    // p-value related values
    public double p_value = 100;
    public int n_random_better_match = -1;
    public int n_total_random_peptides = -1;


    /**
     * Peptide sequence
     */
    public String peptideSequence;

    public Peptide peptide = null;

    public JPeptide jPeptide = null;



    public String spectrum_title;


    /**
     * save additional information as string
     */
    public String eInfo;
}
