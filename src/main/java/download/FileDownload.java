package main.java.download;
import com.alibaba.fastjson.JSON;
import com.alibaba.fastjson.JSONArray;
import com.alibaba.fastjson.JSONObject;
import main.java.raw.RawDataConvert;
import main.java.util.Cloger;
import org.apache.commons.cli.*;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.net.ftp.FTPClient;
import org.apache.commons.net.ftp.FTPFile;
import java.io.*;
import java.lang.reflect.Array;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class FileDownload {

    public FileDownload(){

    }

    public static RawDataConvert rawDataConvert = new RawDataConvert();
    public static boolean do_raw_data_conversion = false;

    private final String proxi_base = "http://proteomecentral.proteomexchange.org/api/proxi/v0.1/datasets/";
    private final String usi_base = "http://proteomecentral.proteomexchange.org/api/proxi/v0.1/spectra?resultType=full&usi=";
    /**
     * Cannot use ftp://massive.ucsd.edu/
     * https://blog.csdn.net/qq_36338455/article/details/86526236
     */
    private final String massive_base = "massive.ucsd.edu";
    private final String jpost_base = "ftp.biosciencedbc.jp";

    public boolean only_show_files = false;
    /**
     * This is useful when process PDC dataset
     */
    public String dataset_id = "";

    public static void main(String[] args) throws IOException, ParseException {

        Cloger.getInstance().logger.info("Start analysis");
        Cloger.getInstance().logger.info(StringUtils.join(args," "));

        Options options = new Options();

        options.addOption("i", true, "A PDC manifest file or a dataset ID");
        options.addOption("o",true,"Output folder");
        options.addOption("f",true,"File pattern to download.");
        options.addOption("c",true,"The number of CPUs. Default is all available CPUs.");
        options.addOption("show",false,"Show all files without downloading any file.");

        options.addOption("h", false, "Help");

        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = parser.parse(options, args);
        if (cmd.hasOption("h") || cmd.hasOption("help") || args.length == 0) {
            HelpFormatter f = new HelpFormatter();
            f.printHelp("Options", options);
            System.exit(0);
        }

        String input = cmd.getOptionValue("i");
        String outdir = "." + File.separator;
        if(cmd.hasOption("o")) {
            outdir = cmd.getOptionValue("o");
        }

        int ncpu = 0;
        if(cmd.hasOption("c")){
            ncpu = Integer.parseInt(cmd.getOptionValue("c"));
        }


        long start_time = System.currentTimeMillis();

        FileDownload fileDownload = new FileDownload();
        if(cmd.hasOption("show")){
            fileDownload.only_show_files = true;
        }

        String file_pattern = "";
        if(cmd.hasOption("f")) {
            file_pattern = cmd.getOptionValue("f");
        }

        fileDownload.download_from_public_db(input,outdir,ncpu,file_pattern);

        long ctime = System.currentTimeMillis();
        double t = 1.0*(ctime  - start_time)/1000.0/60.0;
        String out = String.format("%.2f",t);
        out = out + " min";
        Cloger.getInstance().logger.info("Total time used:" + out);


    }

    /**
     * Download data from public databases, such as PDC, PRIDE and MassIVE
     * @param dataset dataset ID
     * @param outdir output folder
     * @param ncpu the number of threads
     * @param file_pattern file type, file suffix, such as raw, mgf, mzML et. al.
     * @return the folder for the downloaded files.
     */
    public String download_from_public_db(String dataset, String outdir, int ncpu, String file_pattern){

        String data_out_dir = outdir;
        File F = new File(dataset);
        if(F.isFile()){
            try {
                data_out_dir = download_from_pdc(dataset,outdir,ncpu);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }else{
            if(dataset.startsWith("PXD")) {
                data_out_dir = download_from_px(dataset, outdir, ncpu, file_pattern);
            }else if(dataset.startsWith("MSV")) {
                data_out_dir = download_from_massive(dataset, outdir, ncpu, file_pattern);
            }else if(dataset.startsWith("JPST")){
                data_out_dir = download_from_jpost(dataset, outdir, ncpu, file_pattern);
            }else{
                Cloger.getInstance().logger.error("Invalid dataset ID!");
                System.exit(1);
            }
        }
        return data_out_dir;
    }


    private String download_from_pdc(String manifest_file_tsv, String out_dir, int n_cpu) throws IOException {
        BufferedReader bReader = new BufferedReader(new FileReader(new File(manifest_file_tsv)));

        String[] head = bReader.readLine().trim().split("\t");
        HashMap<String,Integer> column_name2index = new HashMap<>();
        for(int i=0;i<head.length;i++){
            column_name2index.put(head[i],i);
        }

        ArrayList<String> urls = new ArrayList<>();
        ArrayList<String> names = new ArrayList<>();

        String study_name = "";
        String study_id = "";

        String line;
        // GB
        double total_data_size = 0;
        while((line = bReader.readLine())!=null){
            line = line.trim();
            String d[] = line.split("\t");

            //
            String file_name = d[column_name2index.get("File Name")];
            String url = d[column_name2index.get("File Download Link")];
            urls.add(url);
            names.add(file_name);

            if(column_name2index.containsKey("Study Name")){
                study_name = d[column_name2index.get("Study Name")];
            }


            if(column_name2index.containsKey("PDC Study ID")){
                study_id = d[column_name2index.get("PDC Study ID")];
            }

            if(column_name2index.containsKey("File Size (in bytes)")){
                total_data_size = total_data_size + Double.parseDouble(d[column_name2index.get("File Size (in bytes)")])/1024.0/1024.0/1024.0;
            }
        }

        bReader.close();

        if(total_data_size>0){
            Cloger.getInstance().logger.info("Total data size:"+String.format("%.2f",total_data_size)+" GB");
        }

        this.dataset_id = study_id;

        Cloger.getInstance().logger.info("Study name: "+study_name);
        Cloger.getInstance().logger.info("Study ID: "+study_id);
        Cloger.getInstance().logger.info("Total files: "+urls.size());

        if(this.only_show_files) {
            int i = 0;
            for(String url:urls){
                i = i + 1;
                System.out.println(i+"\t"+url);
            }

        }else{
            if(!study_id.isEmpty()){
                out_dir = out_dir + File.separator + study_id;
            }

            File F = new File(out_dir);
            if(F.exists()){
                F.mkdirs();
            }
            download(urls,names,out_dir,n_cpu);
        }

        return out_dir;

    }

    private String download_from_px(String dataset_id, String out_dir, int n_cpu, String file_pattern){

        this.dataset_id = dataset_id;

        ArrayList<String> urls = new ArrayList<>();
        ArrayList<String> names = new ArrayList<>();

        String json = "";

        try {
            json = IOUtils.toString(new URL(this.proxi_base+dataset_id), StandardCharsets.UTF_8);
        } catch (IOException e) {
            Cloger.getInstance().logger.error("DataSet ID is not valid or check your internet connection!");
            e.printStackTrace();
            System.exit(1);
        }

        if (json.contains("datasetFiles")) {
            JSONObject jsonObject = (JSONObject) JSONArray.parse(json);
            for (Map.Entry<String, Object> entry : jsonObject.entrySet()) {
                //System.out.println(entry.getKey());
                if (entry.getKey().equalsIgnoreCase("datasetFiles")) {
                    List<px_data_file> fs = JSON.parseArray(entry.getValue().toString(), FileDownload.px_data_file.class);
                    for (FileDownload.px_data_file f: fs) {
                        if (f.getValue() != null && !f.getValue().isEmpty()) {
                            String url = f.getValue();
                            if (file_pattern.isEmpty()) {
                                urls.add(url);
                            } else {
                                // a pattern for selecting a set of datasets
                                Pattern pattern = Pattern.compile(file_pattern);
                                Matcher matcher = pattern.matcher(url);
                                if(matcher.find()){
                                    urls.add(url);
                                }
                            }
                        }
                    }
                }

            }
        }else{
            Cloger.getInstance().logger.error("There is no file found in this dataset:"+dataset_id);
            System.exit(1);
        }

        if(urls.size()<=0){
            Cloger.getInstance().logger.error("There is no file found in this dataset:"+dataset_id);
            System.exit(1);
        }

        Cloger.getInstance().logger.info("Dataset ID: "+dataset_id);
        Cloger.getInstance().logger.info("Total files: "+urls.size());

        if(this.only_show_files) {
            int i = 0;
            for(String url:urls){
                i = i + 1;
                System.out.println(i+"\t"+url);
            }

        }else{
            out_dir = out_dir + File.separator + dataset_id;

            File F = new File(out_dir);
            if (F.exists()) {
                F.mkdirs();
            }
            download(urls, names, out_dir, n_cpu);
        }
        return out_dir;

    }

    private String download_from_massive(String dataset_id, String out_dir, int n_cpu, String file_pattern){
        String new_out_dir = download_from_ftp(dataset_id,"", out_dir,n_cpu,file_pattern,this.massive_base);
        return new_out_dir;
    }

    private String download_from_jpost(String dataset_id, String out_dir, int n_cpu, String file_pattern){
        String new_out_dir = download_from_ftp(dataset_id,"archive/jpostrepos",out_dir,n_cpu,file_pattern,this.jpost_base);
        return new_out_dir;
    }


    public String download_from_ftp(String dataset_id, String host_extra, String out_dir, int n_cpu, String file_pattern, String ftp_base){

        this.dataset_id = dataset_id;

        ArrayList<String> urls = new ArrayList<>();
        ArrayList<String> names = new ArrayList<>();

        String remoteDirectory;
        if(host_extra.isEmpty()) {
            remoteDirectory = "/" + dataset_id; // change if different
        }else{
            remoteDirectory = "/" + host_extra + "/" + dataset_id; // change if different
        }

        double total_data_size = get_file_list_from_ftp(ftp_base, remoteDirectory, file_pattern, urls);


        if(urls.size()<=0){
            Cloger.getInstance().logger.error("There is no file found in this dataset:"+dataset_id);
            System.exit(1);
        }

        Cloger.getInstance().logger.info("Dataset ID: "+dataset_id);
        Cloger.getInstance().logger.info("Total files: "+urls.size());
        Cloger.getInstance().logger.info("Total data size: "+String.format("%.2fGB",total_data_size));


        if(this.only_show_files) {
            int i = 0;
            for(String url:urls){
                i = i + 1;
                System.out.println(i+"\t"+url);
            }

        }else{
            out_dir = out_dir + File.separator + dataset_id;

            File F = new File(out_dir);
            if (F.exists()) {
                F.mkdirs();
            }
            download(urls, names, out_dir, n_cpu);
        }

        return out_dir;

    }

    public static double get_file_list_from_ftp(String ftp_base, String remoteDirectory, String file_pattern, ArrayList<String> urls){

        ArrayList<String> base_ftp_file_names = new ArrayList<>();
        ArrayList<Double> data_file_sizes = new ArrayList<>();

        get_file_list_from_ftp(ftp_base, remoteDirectory, base_ftp_file_names, data_file_sizes);

        double total_data_size = 0.0;
        if(base_ftp_file_names.size()>0){
            for(int i=0;i<base_ftp_file_names.size();i++){
                String name = base_ftp_file_names.get(i);

                if (file_pattern.isEmpty()) {
                    urls.add("ftp://"+ftp_base+"/"+name);
                    total_data_size = total_data_size + data_file_sizes.get(i);
                } else {
                    // a pattern for selecting a set of datasets
                    Pattern pattern = Pattern.compile(file_pattern);
                    Matcher matcher = pattern.matcher(name);
                    if(matcher.find()){
                        urls.add("ftp://"+ftp_base+"/"+name);
                        total_data_size = total_data_size + data_file_sizes.get(i);
                    }
                }

            }
        }
        return total_data_size;
    }

    public static void get_file_list_from_ftp(String ftp_base, String remoteDirectory,
                                              ArrayList<String> base_ftp_file_names,
                                              ArrayList<Double> data_file_sizes){

        FTPClient client = new FTPClient();
        String host = ftp_base;
        Integer ftpPort = 21; // default ftp port, change your port number

        try {
            client.connect(host);
        } catch (IOException e) {
            e.printStackTrace();
        }
        client.enterLocalPassiveMode();
        try {
            client.login("anonymous", "");
        } catch (IOException e) {
            e.printStackTrace();
        }

        try {
            list_all_files(client,remoteDirectory,"",base_ftp_file_names,data_file_sizes);
        } catch (IOException e) {
            e.printStackTrace();
        }

        // log out and disconnect from the server
        if(client.isConnected()){
            try {
                client.logout();
            } catch (IOException e) {
                e.printStackTrace();
            }
            try {
                client.disconnect();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    private static void list_all_files(FTPClient ftp_client, String parent_dir, String current_dir, ArrayList<String> target_files, ArrayList<Double> target_file_sizes) throws IOException {
        String dir_to_list = parent_dir;
        if(!current_dir.equals("")){
            dir_to_list = dir_to_list + "/" + current_dir;
        }
        FTPFile[] subFiles = ftp_client.listFiles(dir_to_list);
        if(subFiles != null && subFiles.length >0){
            for(FTPFile aFile : subFiles){
                String current_file_name = aFile.getName();
                if(current_file_name.equals(".") || current_file_name.equals("..")){
                    // skip parent directory and current directory
                    continue;
                }else{
                    if(aFile.isDirectory()){
                        list_all_files(ftp_client,dir_to_list,current_file_name,target_files,target_file_sizes);
                    }else{
                        //System.out.println(parent_dir+"\t"+current_dir+"\t"+current_file_name+"\t"+aFile.getSize());
                        String file_path;
                        if(current_dir.isEmpty()){
                            file_path = parent_dir +"/" + current_file_name;
                        }else{
                            file_path = parent_dir +"/" + current_dir +"/" + current_file_name;
                        }
                        target_files.add(file_path);
                        // Unit in GB
                        target_file_sizes.add(aFile.getSize()/1024.0/1024.0/1024.0);
                    }
                }
            }
        }
    }

    private void download(ArrayList<String> urls, ArrayList<String> names, String out_dir, int n_cpu){

        Cloger.getInstance().logger.info("Downloading "+urls.size()+" files");

        boolean contain_raw_ms_file = false;
        for(String url: urls){
            if(url.toLowerCase().endsWith(".raw")){
                contain_raw_ms_file = true;
                break;
            }
        }

        if(do_raw_data_conversion && contain_raw_ms_file && rawDataConvert.convert_bin.isEmpty()){
            Cloger.getInstance().logger.error("The path of ThermoRawFileParser.exe is not provided!");
            System.exit(1);
        }

        // if(do_raw_data_conversion &&
        //        contain_raw_ms_file &&
        //        rawDataConvert.convert_bin.isEmpty() &&
        //        rawDataConvert.convert_method.equalsIgnoreCase("ThermoRawFileParser")){

        //    File dirF = new File(out_dir);
        //    String tool_dir = dirF.getParent();
        //    if(tool_dir.equalsIgnoreCase(".")){
        //        tool_dir = out_dir;
        //    }
        //    String tool_bin = rawDataConvert.download_ThermoRawFileParser(tool_dir);
        //    rawDataConvert.convert_bin = tool_bin;
        // }

        if (n_cpu == 0) {
            n_cpu = Runtime.getRuntime().availableProcessors();
        }

        if(n_cpu >= 10){
            n_cpu = 10;
        }

        Cloger.getInstance().logger.info("Used CPUs: " + n_cpu);

        ExecutorService fixedThreadPool = Executors.newFixedThreadPool(n_cpu);

        for(int i=0;i<urls.size();i++){
            String out_file_path = "";
            if(names.size()<=0){
                File F = new File(urls.get(i));
                out_file_path = out_dir + File.separator + F.getName();
            }else{
                out_file_path = out_dir + File.separator + names.get(i);
            }
            fixedThreadPool.execute(new DownloadWorker(urls.get(i),out_file_path));
        }

        fixedThreadPool.shutdown();

        try {
            fixedThreadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        Cloger.getInstance().logger.info("File downloading finished");
    }

    public static void download_single_file(String url, String out_file){
        Cloger.getInstance().logger.info("Download file "+url+ " to "+out_file);
        try {
            FileUtils.copyURLToFile(new URL(url),new File(out_file));
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
        }

        Cloger.getInstance().logger.info("Download file "+url+ " to "+out_file+" done");

        if(do_raw_data_conversion){
            if(out_file.toLowerCase().endsWith(".raw")) {
                File F = new File(out_file);
                if(F.exists()){
                    String out_dir = F.getParent();
                    if(out_dir.isEmpty()){
                        out_dir = "." + File.separator;
                    }
                    rawDataConvert.convert(out_file, out_dir);
                }
            }
        }


    }

    public static class px_data_file{

        public String value;
        public String getValue(){
            return value;
        }
        public void setValue(String value){
            this.value = value;
        }
    }

    public String download_usi(String usi_url, String out_dir) throws IOException {

        String json = "";
        String raw_json = "";

        try {
            raw_json = IOUtils.toString(new URL(this.usi_base+usi_url), StandardCharsets.UTF_8);
            json = format_json(raw_json);
        } catch (IOException e) {
            Cloger.getInstance().logger.error("USI is not valid");
            e.printStackTrace();
            System.exit(1);
        }
        ArrayList<Double> intList = new ArrayList<>();
        ArrayList<Double> mzList = new ArrayList<>();
        double precursor_mz = 0;
        int precursor_charge = 0;
        boolean has_precursor_charge = false;
        int scan = 0;
        boolean has_scan = false;
        String spectrum_name = "";
        boolean has_spectrum_name = false;
        String unmodified_peptide_sequence = "";
        boolean has_unmodified_peptide_sequence = false;

        if (json.contains("intensities")) {
            JSONObject jsonObject = (JSONObject) JSONArray.parse(json);
            for (Map.Entry<String, Object> entry : jsonObject.entrySet()) {
                if (entry.getKey().equalsIgnoreCase("intensities")) {
                    intList = (ArrayList<Double>) JSON.parseArray(entry.getValue().toString(), Double.class);
                }else if (entry.getKey().equalsIgnoreCase("mzs")) {
                    mzList = (ArrayList<Double>) JSON.parseArray(entry.getValue().toString(), Double.class);
                }else if(entry.getKey().equalsIgnoreCase("attributes")){
                    List<usi_json_attributes> fs = JSON.parseArray(entry.getValue().toString(), FileDownload.usi_json_attributes.class);
                    for (FileDownload.usi_json_attributes f: fs) {
                        if(f.accession.equalsIgnoreCase("MS:1000827")){
                            // "name": "isolation window target m/z"
                            precursor_mz = Double.parseDouble(f.value);
                        }else if(f.accession.equalsIgnoreCase("MS:1000041")){
                            // "name": "charge state"
                            precursor_charge = Integer.parseInt(f.value);
                            has_precursor_charge = true;
                        }else if(f.accession.equalsIgnoreCase("MS:1008025")){
                            // "name": "scan number"
                            scan = Integer.parseInt(f.value);
                            has_scan = true;
                        }else if(f.accession.equalsIgnoreCase("MS:1003061")){
                            // "name": "spectrum name"
                            spectrum_name = f.value;
                            has_spectrum_name = true;
                        }else if(f.accession.equalsIgnoreCase("MS:1000888")){
                            has_unmodified_peptide_sequence = true;
                            unmodified_peptide_sequence = f.value;
                        }
                    }
                }
            }
        }else{
            Cloger.getInstance().logger.error("Invalid USI json format:"+json);
            System.exit(1);
        }

        StringBuilder mgf = new StringBuilder();
        if(mzList.size()>=1 && mzList.size() == intList.size()){

            mgf.append("BEGIN IONS\n");
            mgf.append("TITLE=").append(usi_url).append("\n");
            mgf.append("PEPMASS=").append(precursor_mz).append("\n");
            if(has_precursor_charge) {
                mgf.append("CHARGE=").append(precursor_charge).append("+").append("\n");
            }else{
                Cloger.getInstance().logger.warn("Precursor charge is not present in the USI.");
            }
            if(has_scan){
                mgf.append("SCANS=").append(scan).append("\n");
            }
            for(int i=0;i<mzList.size();i++){
                mgf.append(mzList.get(i)).append(" ").append(intList.get(i)).append("\n");
            }
            mgf.append("END IONS\n");

        }else{
            Cloger.getInstance().logger.error("Invalid USI format!");
            System.exit(1);
        }

        String usi_dir = out_dir + File.separator + "usi";
        File F = new File(usi_dir);
        if(!F.isDirectory()){
            F.mkdirs();
        }
        String out_json_file = usi_dir + File.separator + "usi.json";
        BufferedWriter jWriter = new BufferedWriter(new FileWriter(out_json_file));
        jWriter.write(raw_json+"\n");
        jWriter.close();

        String out_mgf_file = usi_dir + File.separator + "usi.mgf";
        BufferedWriter mWriter = new BufferedWriter(new FileWriter(out_mgf_file));
        mWriter.write(mgf+"\n");
        mWriter.close();

        return out_mgf_file;
    }

    public static class usi_json_attributes{

        public String value;
        public String name;
        public String accession;
        public String getValue(){
            return value;
        }
        public void setValue(String value){
            this.value = value;
        }

        public String getName() {
            return name;
        }

        public void setName(String name) {
            this.name = name;
        }

        public String getAccession() {
            return accession;
        }

        public void setAccession(String accession) {
            this.accession = accession;
        }
    }

    public String format_json(String json){
        if(json.startsWith("[")){
            json = json.replaceAll("^\\[*","");
            json = json.replaceAll("]*$","");
        }
        return json;
    }

    public static boolean isValidURL(String url) {
        /* Try creating a valid URL */
        try {
            IOUtils.toString(new URL(url), StandardCharsets.UTF_8);
            return true;
        }

        // If there was an Exception
        // while creating URL object
        catch (Exception e) {
            return false;
        }
    }

}
