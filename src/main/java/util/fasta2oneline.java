package main.java.util;

import net.sf.jfasta.FASTAElement;
import net.sf.jfasta.FASTAFileReader;
import net.sf.jfasta.impl.FASTAElementIterator;
import net.sf.jfasta.impl.FASTAFileReaderImpl;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;


public class fasta2oneline {

    public static void main(String[] args) throws IOException {

        String db = args[0];
        String newdb = args[1];
        File dbFile = new File(db);
        FASTAFileReader reader = new FASTAFileReaderImpl(dbFile);
        FASTAElementIterator it = reader.getIterator();

        BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(new File(newdb)));

        while (it.hasNext()) {
            FASTAElement el = it.next();
            el.setLineLength(1);

            bufferedWriter.write(">"+el.getHeader()+"\n"+el.getSequence()+"\n");

        }

        reader.close();
        bufferedWriter.close();


    }

}
