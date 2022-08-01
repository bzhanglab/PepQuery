package main.java.raw;
import main.java.util.Cloger;
import org.apache.commons.cli.*;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.SystemUtils;

import java.io.*;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

public class RawDataConvert {

    /**
     * The tool used for data conversion: ThermoRawFileParser or MSConvert
     * ThermoRawFileParser: <a href="https://github.com/compomics/ThermoRawFileParser">ThermoRawFileParser</a>
     */
    public String convert_method = "ThermoRawFileParser";
    public String convert_bin = "";

    public boolean delete_raw_file_after_conversion = false;

    /**
     * The format of MS/MS data: mzML, mgf
     */
    public String ms_format = "mzML";
    public boolean do_gz = true;

    public static void main(String[] args) throws IOException, ParseException {

        Cloger.getInstance().logger.info("Start analysis");
        Cloger.getInstance().logger.info(StringUtils.join(args," "));

        Options options = new Options();

        options.addOption("i", true, "A raw MS/MS file or a folder");
        options.addOption("o",true,"Output folder");
        options.addOption("f",true,"File format: mgf or mzML");
        options.addOption("c",true,"The number of CPUs. Default is all available CPUs");
        options.addOption("p",true,"Tool path for Raw MS/MS data conversion: for example, /home/test/ThermoRawFileParser/ThermoRawFileParser.exe");
        options.addOption("m",true,"Method for Raw MS/MS data conversion: ThermoRawFileParser");
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

        RawDataConvert rawDataConvert = new RawDataConvert();
        if(cmd.hasOption("f")) {
            rawDataConvert.ms_format = cmd.getOptionValue("f");
        }

        if(cmd.hasOption("m")){
            rawDataConvert.convert_method = cmd.getOptionValue("m");
        }

        rawDataConvert.convert_bin = cmd.getOptionValue("p");

        rawDataConvert.convert(input, outdir);

        long ctime = System.currentTimeMillis();
        double t = 1.0*(ctime  - start_time)/1000.0/60.0;
        String out = String.format("%.2f",t);
        out = out + " min";
        Cloger.getInstance().logger.info("Total time used:" + out);

    }

    public boolean convert(String raw_file, String out_dir){
        boolean convert_finished_without_error = false;
        String cmd = "";

        File rawF = new File(raw_file);
        String file_name = rawF.getName();
        if(file_name.endsWith("raw")){
            file_name = file_name.replaceAll("\\.raw$","");
        }else if(file_name.endsWith("RAW")){
            file_name = file_name.replaceAll("\\.RAW$","");
        }else{
            file_name = file_name.replaceAll("\\.[^.]+$","");
        }


        if(this.convert_method.equalsIgnoreCase("ThermoRawFileParser")){
            if(SystemUtils.IS_OS_LINUX){
                cmd = "mono "+this.convert_bin;
            }else if(SystemUtils.IS_OS_WINDOWS){
                cmd = this.convert_bin;
            }else {
                Cloger.getInstance().logger.error("OS is not supported yet!");
                System.exit(1);
            }

            if(this.do_gz){
                cmd = cmd + " -g -i "+ raw_file;
            }else{
                cmd = cmd + " -i "+ raw_file;
            }

            if(this.ms_format.equalsIgnoreCase("mzML")){
                // 2 for indexed mzML
                cmd = cmd + " -f 2 -L 1- -b " + out_dir+File.separator+file_name+".mzML";
            }else if(this.ms_format.equalsIgnoreCase("mgf")){
                // 0 for MGF
                cmd = cmd + " -f 0 -L 2 -b "+out_dir+File.separator+file_name+".mgf";
            }else{
                Cloger.getInstance().logger.error("MS/MS format is not supported yet:"+this.ms_format);
                System.exit(1);
            }
        }else if(this.convert_method.equalsIgnoreCase("MSConvert")){

        }else{
            Cloger.getInstance().logger.error("Method is not supported yet:"+this.convert_method);
            System.exit(1);
        }

        Cloger.getInstance().logger.info("Convert "+raw_file + ": "+ cmd);
        boolean conversion_pass = run_cmd(cmd);
        if(conversion_pass && this.delete_raw_file_after_conversion){
            if(!rawF.delete()){
                Cloger.getInstance().logger.warn("File deletion failed:"+raw_file);
            }
        }

        return convert_finished_without_error;
    }

    private boolean run_cmd(String cmd){
        boolean conversion_pass = true;
        Runtime rt = Runtime.getRuntime();
        Process p;
        try {
            p = rt.exec(cmd);
        } catch (IOException e) {
            conversion_pass = false;
            throw new RuntimeException(e);
        }
        StreamLog errorLog = new StreamLog(p.getErrorStream(), this.convert_method+" => Error:", true);
        StreamLog stdLog = new StreamLog(p.getInputStream(), this.convert_method+" => Message:", true);

        errorLog.start();
        stdLog.start();

        try {
            int exitValue = p.waitFor();
            if (exitValue != 0) {
                conversion_pass = false;
                Cloger.getInstance().logger.error("MS data conversion error:" + exitValue);
            }
        } catch (InterruptedException e) {
            conversion_pass = false;
            throw new RuntimeException(e);
        }

        try {
            errorLog.join();
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        }
        try {
            stdLog.join();
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        }

        return conversion_pass;
    }

    /**
    public String download_ThermoRawFileParser(String out_dir){
        Cloger.getInstance().logger.info("Download ThermoRawFileParser for raw MS/MS data conversion");
        String link = "https://github.com/compomics/ThermoRawFileParser/releases/download/v1.3.4/ThermoRawFileParser.zip";
        String out_file = out_dir + File.separator + "ThermoRawFileParser.zip";
        FileDownload.download_single_file(link,out_file);
        try {
            unzip(out_file,out_dir+File.separator+"ThermoRawFileParser");
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        String tool_bin = out_dir + File.separator + "ThermoRawFileParser" + File.separator + "ThermoRawFileParser.exe";
        this.convert_bin = tool_bin;
        Cloger.getInstance().logger.info("ThermoRawFileParser path:"+tool_bin);
        return tool_bin;
    }
     **/

    public void unzip(String zip_file, String out_dir) throws IOException {
        // reference: https://www.baeldung.com/java-compress-and-uncompress

        File dirF = new File(out_dir);
        if(!dirF.exists()){
            dirF.mkdirs();
        }
        byte[] buffer = new byte[1024];
        ZipInputStream zis = new ZipInputStream(new FileInputStream(zip_file));
        ZipEntry zipEntry = zis.getNextEntry();

        while (zipEntry != null) {
            File newFile = new File(dirF, zipEntry.getName());

            String outDirPath = dirF.getCanonicalPath();
            String outFilePath = newFile.getCanonicalPath();

            if (!outFilePath.startsWith(outDirPath + File.separator)) {
                throw new IOException("Entry is outside of the target dir: " + zipEntry.getName());
            }
            if (zipEntry.isDirectory()) {
                if (!newFile.isDirectory() && !newFile.mkdirs()) {
                    throw new IOException("Failed to create directory " + newFile);
                }
            } else {
                File parent = newFile.getParentFile();
                if (!parent.isDirectory() && !parent.mkdirs()) {
                    throw new IOException("Failed to create directory " + parent);
                }

                FileOutputStream fos = new FileOutputStream(newFile);
                int len;
                while ((len = zis.read(buffer)) > 0) {
                    fos.write(buffer, 0, len);
                }
                fos.close();
            }
            zipEntry = zis.getNextEntry();
        }
        zis.closeEntry();
        zis.close();
    }




}
