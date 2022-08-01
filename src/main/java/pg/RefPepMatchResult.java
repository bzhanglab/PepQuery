package main.java.pg;

import java.util.ArrayList;

public class RefPepMatchResult {
    String spectrum_title;
    int n_matched_peptides = 0;
    int n_matched_peptides_with_better_score = 0;
    double score = -1;
    ArrayList<String> out_data = new ArrayList<>();

    ArrayList<String> psm_annotation_data = new ArrayList<>();

}
