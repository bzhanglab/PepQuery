package main.java.raw;

import java.io.BufferedReader;

import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class StreamLog extends Thread {
    InputStream is;
    String type;
    boolean isOut;

    public StreamLog(InputStream is, String type, boolean isOut) {
        this.is = is;
        this.type = type;
        this.isOut = isOut;
    }

    public void run() {
        try {
            InputStreamReader isr = new InputStreamReader(is);
            BufferedReader br = new BufferedReader(isr);
            String line = null;
            while ((line = br.readLine()) != null) {
                if (this.isOut) {
                    System.out.println(type + ">" + line);
                }
            }
        } catch (IOException ioe) {
            ioe.printStackTrace();
            System.exit(1);
        }
    }
}
