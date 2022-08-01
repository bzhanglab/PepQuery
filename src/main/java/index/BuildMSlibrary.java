package main.java.index;

import main.java.download.FileDownload;
import main.java.util.Cloger;
import org.apache.commons.cli.*;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.zip.GZIPOutputStream;

/**
 * Build MS/MS library for fast peptide matching.
 */
public class BuildMSlibrary {
    private final int fold = 10;
    public int ncpu = 0;

    /**
     * If there is an error when compress a file, this parameter will need to set as true.
     */
    public static boolean compress_file_error = false;

    private static Logger logger = LogManager.getLogger(BuildMSlibrary.class.getName());


    public static void main(String[] args) throws IOException, ParseException {

        logger.info("Start analysis");
        logger.info(StringUtils.join(args," "));

        Options options = new Options();


        options.addOption("i", true, "MS/MS data file or folder");
        options.addOption("d", true, "Dataset ID");
        options.addOption("o",true,"Output folder");
        options.addOption("c",true,"The number of CPUs. Default is all available CPUs.");
        options.addOption("f",true,"File type (file suffix) to download.");
        options.addOption("p",true,"Tool path for Raw MS/MS data conversion: for example, /home/test/ThermoRawFileParser/ThermoRawFileParser.exe");
        options.addOption("m",true,"Method for Raw MS/MS data conversion: ThermoRawFileParser");
        options.addOption("r",false,"Delete raw file after data conversion");
        options.addOption("show",false,"Show all files without downloading any file.");

        options.addOption("h", false, "Help");

        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = parser.parse(options, args);
        if (cmd.hasOption("h") || cmd.hasOption("help") || args.length == 0) {
            HelpFormatter f = new HelpFormatter();
            // System.out.println(version);
            f.printHelp("Options", options);
            System.exit(0);
        }


        String outdir = cmd.getOptionValue("o");

        int ncpu = 0;
        if(cmd.hasOption("c")){
            ncpu = Integer.parseInt(cmd.getOptionValue("c"));
        }

        long start_time = System.currentTimeMillis();

        BuildMSlibrary buildMSlibrary = new BuildMSlibrary();
        buildMSlibrary.ncpu = ncpu;


        FileDownload.do_raw_data_conversion = true;
        if(cmd.hasOption("m")){
            FileDownload.rawDataConvert.convert_method = cmd.getOptionValue("m");
        }

        if(cmd.hasOption("p")){
            FileDownload.rawDataConvert.convert_bin = cmd.getOptionValue("p");
        }

        if(cmd.hasOption("r")){
            FileDownload.rawDataConvert.delete_raw_file_after_conversion = true;
        }

        FileDownload.do_raw_data_conversion = true;
        FileDownload.rawDataConvert.delete_raw_file_after_conversion = true;
        FileDownload.rawDataConvert.ms_format = "mgf";
        FileDownload.rawDataConvert.do_gz = false;

        if(cmd.hasOption("i")){
            String msfile = cmd.getOptionValue("i");
            buildMSlibrary.buildLibrary(msfile, outdir);
        }else{
            if(cmd.hasOption("d")) {
                String dataset_id = cmd.getOptionValue("d");
                String file_type = "";
                if (cmd.hasOption("f")) {
                    file_type = cmd.getOptionValue("f");
                }
                boolean only_show_data_file = false;
                if (cmd.hasOption("show")) {
                    only_show_data_file = true;
                }
                buildMSlibrary.buildLibrary(dataset_id, outdir, file_type, only_show_data_file);
            }
        }


        long ctime = System.currentTimeMillis();
        double t = 1.0*(ctime  - start_time)/1000.0/60.0;
        String out = String.format("%.2f",t);
        out = out + " min";
        Cloger.getInstance().logger.info("Total time used:" + out);


    }


    /**
     * Build MS/MS index from public dataset.
     * @param dataset_id
     * @param outdir
     * @param file_type
     * @param only_show_data_file
     */
    public String buildLibrary(String dataset_id, String outdir, String file_type, boolean only_show_data_file){
        FileDownload fileDownload = new FileDownload();
        fileDownload.only_show_files = only_show_data_file;
        String data_dir = outdir + File.separator + "download";
        String ms_file = fileDownload.download_from_public_db(dataset_id,data_dir,ncpu,file_type);
        String index_dir = outdir + File.separator + fileDownload.dataset_id +File.separator + "index";

        if(!only_show_data_file) {
            try {
                buildLibrary(ms_file, index_dir);
            } catch (IOException e) {
                e.printStackTrace();
            }
            Cloger.getInstance().logger.info("Index directory for dataset " + fileDownload.dataset_id + ": " + index_dir);
        }
        return index_dir;
    }


    public String buildLibrary(String msfile, String outdir) throws IOException {
        File msFile = new File(msfile);
        ArrayList<File> spectraList = new ArrayList<>();
        if(msFile.isDirectory()){
            //File mslist[] = msFile.listFiles();
            List<File> mslist = new ArrayList<>();
            try (Stream<Path> walk = Files.walk(Paths.get(msfile))) {

                mslist = walk.filter(Files::isRegularFile).map(path -> path.toFile()).collect(Collectors.toList());


            } catch (IOException e) {
                e.printStackTrace();
            }
            for(int i=0;i<mslist.size();i++){

                if(mslist.get(i).getName().toLowerCase().endsWith(".mzml") ||
                        mslist.get(i).getName().toLowerCase().endsWith(".mzml.gz") ||
                        mslist.get(i).getName().toLowerCase().endsWith(".mgf") ||
                        mslist.get(i).getName().toLowerCase().endsWith(".mgf.gz") ||
                        mslist.get(i).getName().toLowerCase().endsWith(".mzxml") ||
                        mslist.get(i).getName().toLowerCase().endsWith(".mzxml.gz")
                ){
                    spectraList.add(mslist.get(i));
                    Cloger.getInstance().logger.info("Add "+mslist.get(i).getAbsolutePath());
                }

            }
        }else{
            spectraList.add(msFile);
        }

        File OD = new File(outdir);
        if(!OD.isDirectory()){
            OD.mkdirs();
        }

        String summary_file = outdir + File.separator + "summary.txt";
        BufferedWriter sReader = new BufferedWriter(new FileWriter(summary_file));
        Date objDate = new Date();
        sReader.write("date="+objDate.toString()+"\n");

        Cloger.getInstance().logger.info("Total MS/MS files: " + spectraList.size());
        sReader.write("total_ms_file="+spectraList.size()+"\n");

        if (this.ncpu == 0) {
            ncpu = Runtime.getRuntime().availableProcessors();
        }
        // reduce IO pressure
        if(ncpu >= 8){
            ncpu = 8;
        }

        Cloger.getInstance().logger.info("Used CPUs: " + ncpu);


        ExecutorService fixedThreadPool = Executors.newFixedThreadPool(ncpu);

        ArrayList<String> filedir = new ArrayList<>();

        String tmp_dir = outdir + File.separator + "tmp";

        for(int k=0;k<spectraList.size();k++) {
            String fullMsFilePath = spectraList.get(k).getAbsolutePath();
            String k_outdir = tmp_dir + File.separator + spectraList.get(k).getName() + "_" + k;
            filedir.add(k_outdir);
            fixedThreadPool.execute(new IndexWorker(fullMsFilePath, k_outdir, this.fold));
            sReader.write("ms_file="+fullMsFilePath+"\n");
        }

        fixedThreadPool.shutdown();

        try {
            fixedThreadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }


        long n_spectra = 0;
        for( String dir : filedir){
            Cloger.getInstance().logger.info("Process "+dir);
            File FD = new File(dir);
            File[] files = FD.listFiles();
            ExecutorService fixedThreadPool2 = Executors.newFixedThreadPool(ncpu);

            for(File f: files){
                String filename = f.getName();
                if(filename.endsWith(".mgf")){
                    fixedThreadPool2.execute(new CombineWorker(f, outdir));
                    /***
                    String file_path = outdir + File.separator + filename;
                    File MF = new File(file_path);
                    if(!MF.isFile()){
                        MF.createNewFile();
                    }
                    FileWriter fw = new FileWriter(MF,true);
                    BufferedWriter bw = new BufferedWriter(fw);

                    BufferedReader br = new BufferedReader(new FileReader(f));
                    String line;
                    while((line=br.readLine())!=null){
                        line = line.trim();
                        if(line.startsWith("BEGIN")){
                            n_spectra++;
                        }
                        bw.write(line+"\n");
                    }
                    br.close();

                    bw.close();***/
                }
                ///f.delete();

            }

            fixedThreadPool2.shutdown();

            try {
                fixedThreadPool2.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }

            String count_file = dir+"/count.txt";
            File Count_File = new File(count_file);
            BufferedReader nReader = new BufferedReader(new FileReader(Count_File));
            String n_msms = nReader.readLine().trim();
            n_spectra = n_spectra + Long.valueOf(n_msms);
            nReader.close();
            Count_File.delete();

            FD.delete();
        }

        File TMP = new File(tmp_dir);
        TMP.delete();
        Cloger.getInstance().logger.info("Total spectra: "+ n_spectra);
        sReader.write("total_msms_spectra="+n_spectra+"\n");
        sReader.close();

        Cloger.getInstance().logger.info("Compress index files ...");
        gzfiles(outdir,".mgf",outdir);
        Cloger.getInstance().logger.info("Compress index files done");

        if(!compress_file_error) {
            Cloger.getInstance().logger.info("Delete all .mgf index files, keep .mgf.gz index files ...");
            File DIR = new File(outdir);
            File[] all_files = DIR.listFiles();
            boolean delete_with_error = false;
            for(File F : all_files){
                if(F.isFile() && F.getName().toLowerCase().endsWith(".mgf")){
                    if(!F.delete()){
                        delete_with_error = true;
                        if(delete_with_error){
                            Cloger.getInstance().logger.error("Error: delete file -> "+F.getAbsolutePath());
                            System.exit(1);
                        }
                    }
                }
            }
            Cloger.getInstance().logger.info("Delete all .mgf index files, keep .mgf.gz index files. Done");

        }
        return outdir;

    }

    /**
     *
     * @param mass
     * @param fold 10 or 100.
     * @return
     */
    public static int getIntMass(double mass, int fold){
        int intMass = (int) Math.round(mass * fold);
        return intMass;
    }

    public void gzfiles(String in_dir, String suffix, String out_dir){
        File DIR = new File(in_dir);
        File[] all_files = DIR.listFiles();

        if (this.ncpu == 0) {
            ncpu = Runtime.getRuntime().availableProcessors();
        }

        Cloger.getInstance().logger.info("Used CPUs: " + ncpu);


        ExecutorService fixedThreadPool = Executors.newFixedThreadPool(ncpu);

        for(File F : all_files){
            if(F.isFile() && F.getName().toLowerCase().endsWith(suffix.toLowerCase())){
                String out_gz_file = out_dir + File.separator + F.getName() + ".gz";
                fixedThreadPool.execute(new CompressFileWorker(F.getAbsolutePath(),out_gz_file));
            }
        }

        fixedThreadPool.shutdown();

        try {
            fixedThreadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        if(compress_file_error){
            Cloger.getInstance().logger.error("Error: file compression!");
            System.exit(1);
        }
    }


    /**
     *
     * @param source
     * @param target
     * @throws IOException
     */
    public static boolean compressGzip(Path source, Path target){

        boolean finished = true;
        try (GZIPOutputStream gos = new GZIPOutputStream(new FileOutputStream(target.toFile()))) {
            Files.copy(source, gos);
        } catch (IOException e) {
            finished = false;
            System.err.println("ERROR: compress file -> " + source.toFile().getAbsolutePath());
            e.printStackTrace();
            System.exit(1);
        }
        return finished;
    }


}
