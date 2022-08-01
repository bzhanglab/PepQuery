package main.java.index;

import java.io.*;
import java.util.Scanner;

public final class CombineWorker implements Runnable{


    private File input_mgf;
    private String out_dir;

    public CombineWorker(File in_mgf, String out_dir){
        this.input_mgf = in_mgf;
        this.out_dir = out_dir;
    }


    @Override
    public void run() {
        String filename = input_mgf.getName();
        if(filename.endsWith(".mgf")){
            String file_path = out_dir + File.separator + filename;
            File MF = new File(file_path);
            if(!MF.isFile()){
                try {
                    MF.createNewFile();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            FileWriter fw = null;
            try {
                fw = new FileWriter(MF,true);
            } catch (IOException e) {
                e.printStackTrace();
            }
            BufferedWriter bw = new BufferedWriter(fw);

            Scanner sc = null;
            try {
                sc = new Scanner(input_mgf);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
            String line;
            while (sc.hasNextLine()) {
                line = sc.nextLine();
                try {
                    bw.write(line+"\n");
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

            try {
                bw.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
            sc.close();
        }
        input_mgf.delete();
    }
}
