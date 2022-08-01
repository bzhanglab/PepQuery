package main.java.pg;

import com.compomics.util.experiment.biology.modifications.Modification;

public class JPTM {

    public int pos = -1;
    public Modification modification = new Modification();
    public JPTM(int position, Modification ptm){
        this.pos = position;
        this.modification = ptm;
    }
}
