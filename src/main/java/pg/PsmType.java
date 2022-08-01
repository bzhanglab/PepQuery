package main.java.pg;

public enum PsmType {

    //if the target input sequence is not peptide, then use this parameter to specify the type of input event. 1=protein,2=DNA,3=VCF,4=BED,5=GTF;
    no_spectrum_match("There is no matched spectrum"),
    spectrum_matched("Find the spectrum in matching"),
    low_score("The match for the PSM is <= CParameter.minScore"),
    high_p_value("P value is high (e.g., > 0.01)"),
    better_ref_pep("n_db >= 1"),
    better_ref_pep_with_mod("n_ptm >= 1"),
    not_top_rank("rank >= 2"),
    confident("confident"),
    invalid_peptide("peptide is not valid, for example, a novel peptide is present in reference database");

    /**
     * The description.
     */
    public final String description;

    /**
     * Constructor.
     *
     * @param description the description of PsmType
     */
    PsmType(String description) {
        this.description = description;
    }
}
