package main.java.util;


import net.sf.jfasta.FASTAElement;
import net.sf.jfasta.FASTAFileReader;
import net.sf.jfasta.impl.FASTAElementIterator;
import net.sf.jfasta.impl.FASTAFileReaderImpl;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * preprocess ensembl protein database: remove proteins with * in sequences
 */
public class preprocessEnsemblDB {

    public static void main(String[] args) throws IOException {

        String db = args[0];
        String newdb = args[1];
        File dbFile = new File(db);
        FASTAFileReader reader = new FASTAFileReaderImpl(dbFile);
        FASTAElementIterator it = reader.getIterator();

        BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(new File(newdb)));

        int n_total = 0;
        int n_remove = 0;
        int n_valid = 0;
        while (it.hasNext()) {
            FASTAElement el = it.next();
            el.setLineLength(1);
            n_total++;
            String acc = el.getHeader().split("\\s+")[0];
            String seq = el.getSequence();
            if(seq.contains("*")){
                n_remove++;
                Cloger.getInstance().logger.info("Remove: "+acc + "\t"+ seq);
            }else{
                n_valid++;
                bufferedWriter.write(">"+acc+"\n"+seq+"\n");
            }



        }

        reader.close();
        bufferedWriter.close();

        Cloger.getInstance().logger.info("Total proteins:"+n_total);
        Cloger.getInstance().logger.info("Removed proteins:"+n_remove);
        Cloger.getInstance().logger.info("Valid proteins:"+n_valid);

    }

}
