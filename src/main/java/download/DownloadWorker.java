package main.java.download;

public class DownloadWorker implements Runnable{


    public String url;
    public String out_file;

    public DownloadWorker(String url_file, String out_file_path){
        this.url = url_file;
        this.out_file = out_file_path;
    }

    @Override
    public void run() {
        FileDownload.download_single_file(url,out_file);
    }
}
