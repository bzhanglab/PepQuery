package main.java.msio;

import main.java.pg.S3Interface;
import main.java.util.Cloger;
import org.apache.commons.io.FileUtils;

import java.io.File;
import java.io.IOException;

public class IndexFileDownloadWorker implements Runnable {


    public String s3_file_path;
    public String out_dir;
    public IndexFileDownloadWorker(String s3_file_path, String out_dir){
        this.s3_file_path = s3_file_path;
        this.out_dir = out_dir;
    }


    @Override
    public void run() {
        download(s3_file_path,out_dir);
    }

    public void download(String objPath, String out_dir){;
        String out_file = out_dir + File.separator + FileUtils.getFile(objPath).getName();
        File file = new File(out_file);
        try {
            file.createNewFile();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        if(!S3Interface.getInstance().download(objPath, out_file)) {
            Cloger.getInstance().logger.warn(objPath + " doesn't exist");
            file.delete();
        }
    }
}
