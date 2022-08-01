package main.java.index;

import java.nio.file.Path;
import java.nio.file.Paths;

public class CompressFileWorker implements Runnable{

    private Path source;
    private Path target;

    public CompressFileWorker(String source_file, String target_file){
        this.source = Paths.get(source_file);
        this.target = Paths.get(target_file);
    }

    @Override
    public void run() {
        boolean finished = BuildMSlibrary.compressGzip(source,target);
        if(!finished){
            BuildMSlibrary.compress_file_error = true;
        }
    }
}
