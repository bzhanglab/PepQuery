package main.java.pg;

public class TargetedPSM {

    public PsmType status = PsmType.no_spectrum_match;
    public String peptide;
    public String spectrum_title;

    public TargetedPSM(){

    }

    public void update_status(PsmType type){
        this.status = type;
    }

    public PsmType get_status(){
        return this.status;
    }

    public boolean is_status(PsmType type){
        if(this.status.equals(type)){
            return true;
        }else{
            return false;
        }
    }

    public static PsmType get_psm_status(int peptide_length, int rank, int n_db, double p_value, int n_ptm){

        boolean p_value_passed = CParameter.passed_p_value_validation(peptide_length, p_value);
        if(rank >= 2){
            return PsmType.not_top_rank;
        }else if(n_db >= 1){
            return PsmType.better_ref_pep;
        }else if(!p_value_passed){
            return PsmType.high_p_value;
        }else if(n_ptm >= 1){
            return PsmType.better_ref_pep_with_mod;
        }else{
            return PsmType.confident;
        }
    }

}



