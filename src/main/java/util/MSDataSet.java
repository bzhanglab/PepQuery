package main.java.util;

import com.alibaba.fastjson.JSONArray;
import com.alibaba.fastjson.JSONObject;
import main.java.pg.S3Interface;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;

import java.io.*;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MSDataSet {


    private String dataset_name = "";
    /**
     * For visualization purpose, may not be unique.
     */
    private String short_name = "";
    private String ms_file = "";
    private String sample_annotation_file = "";
    private String protein_reference_database = "";
    private String pga_annotation_dir = "";
    private String parameter = "";
    private String species = "";
    private String data_type = "";
    private String data_link = "";


    public CParameterSet parameterSet = new CParameterSet();

    public MSDataSet(){

    }

    public String getDataset_name() {
        return dataset_name;
    }

    public void setDataset_name(String dataset_name) {
        this.dataset_name = dataset_name;
    }

    public String getMs_file() {
        return ms_file;
    }

    public void setMs_file(String ms_file) {
        this.ms_file = ms_file;
    }

    public String getSample_annotation_file() {
        return sample_annotation_file;
    }

    public void setSample_annotation_file(String sample_annotation_file) {
        this.sample_annotation_file = sample_annotation_file;
    }

    public String getProtein_reference_database() {
        return protein_reference_database;
    }

    public void setProtein_reference_database(String protein_reference_database) {
        this.protein_reference_database = protein_reference_database;
    }

    public String getPga_annotation_dir() {
        return pga_annotation_dir;
    }

    public void setPga_annotation_dir(String pga_annotation_dir) {
        this.pga_annotation_dir = pga_annotation_dir;
    }

    public String getParameter() {
        return parameter;
    }

    public void setParameter(String parameter) {
        this.parameter = parameter;
    }

    public String getSpecies() {
        return species;
    }

    public void setSpecies(String species) {
        this.species = species;
    }

    public String getData_type() {
        return data_type;
    }

    public void setData_type(String data_type) {
        this.data_type = data_type;
    }

    public String getData_link() {
        return data_link;
    }

    public void setData_link(String data_link) {
        this.data_link = data_link;
    }

    public String getShort_name() {
        return short_name;
    }

    public void setShort_name(String short_name) {
        this.short_name = short_name;
    }

    public static HashMap<String, MSDataSet> load_default_datasets(){

        HashMap<String,CParameterSet> parameterSetHashMap = CParameterSet.load_default_parameter_sets();
        // java.io.InputStream
        InputStream inputStream = MSDataSet.class.getResourceAsStream("/main/resources/msms.json");
        String json = null;
        try {
            json = IOUtils.toString(inputStream, StandardCharsets.UTF_8.name());
        } catch (IOException e) {
            e.printStackTrace();
        }
        //System.out.println(JSON.parse(json).toString());
        JSONObject jsonObject = (JSONObject) JSONArray.parse(json);
        //System.out.println(jsonObject.size());
        HashMap<String,MSDataSet> msDataSetHashMap = new HashMap<>();
        for (Map.Entry<String, Object> entry : jsonObject.entrySet()) {
            //System.out.println(entry.getKey());
            String dataset_id = entry.getKey();
            MSDataSet msDataSet = ((JSONObject) entry.getValue()).toJavaObject(MSDataSet.class);
            if(parameterSetHashMap.containsKey(msDataSet.getParameter())){
                msDataSet.parameterSet = parameterSetHashMap.get(msDataSet.getParameter());
                msDataSetHashMap.put(dataset_id,msDataSet);
            }else{
                Cloger.getInstance().logger.error("Parameter set doesn't exist:"+msDataSet.getParameter());
                System.exit(1);
            }
        }
        Cloger.getInstance().logger.info("Load "+msDataSetHashMap.size()+" MS/MS datasets.");
        return msDataSetHashMap;
    }

    public static void print_default_datasets(boolean show_full_information){
        HashMap<String, MSDataSet> dataSetHashMap = load_default_datasets();
        ArrayList<String> names = new ArrayList<>(dataSetHashMap.keySet());
        Collections.sort(names);

        int i=0;

        HashMap<String,Integer> column2str_length = new HashMap<>();
        HashMap<String,Long> dataset2n_spectra = new HashMap<>();
        HashMap<String,Long> dataset2n_ms_file = new HashMap<>();
        for(String name: names){
            if(column2str_length.containsKey("dataset_name")){
                if(name.length() > column2str_length.get("dataset_name")){
                    column2str_length.put("dataset_name",name.length());
                }
            }else{
                column2str_length.put("dataset_name",name.length());
            }

            if(column2str_length.containsKey("parameter_set")){
                if(dataSetHashMap.get(name).getParameter().length() > column2str_length.get("parameter_set")){
                    column2str_length.put("parameter_set",dataSetHashMap.get(name).getParameter().length());
                }
            }else{
                column2str_length.put("parameter_set",dataSetHashMap.get(name).getParameter().length());
            }

            if(column2str_length.containsKey("species")){
                if(dataSetHashMap.get(name).getSpecies().length() > column2str_length.get("species")){
                    column2str_length.put("species",dataSetHashMap.get(name).getSpecies().length());
                }
            }else{
                column2str_length.put("species",dataSetHashMap.get(name).getSpecies().length());
            }

            if(column2str_length.containsKey("data_type")){
                if(dataSetHashMap.get(name).getData_type().length() > column2str_length.get("data_type")){
                    column2str_length.put("data_type",dataSetHashMap.get(name).getData_type().length());
                }
            }else{
                column2str_length.put("data_type",dataSetHashMap.get(name).getData_type().length());
            }

            if(column2str_length.containsKey("data_link")){
                if(dataSetHashMap.get(name).getData_link().length() > column2str_length.get("data_link")){
                    column2str_length.put("data_link",dataSetHashMap.get(name).getData_link().length());
                }
            }else{
                column2str_length.put("data_link",dataSetHashMap.get(name).getData_link().length());
            }

            if(show_full_information){
                dataset2n_spectra.put(name,get_the_number_of_spectra(dataSetHashMap.get(name)));
                dataset2n_ms_file.put(name,get_the_number_of_ms_files(dataSetHashMap.get(name)));

                if(column2str_length.containsKey("short_name")){
                    if(dataSetHashMap.get(name).getShort_name().length() > column2str_length.get("short_name")){
                        column2str_length.put("short_name",dataSetHashMap.get(name).getShort_name().length());
                    }
                }else{
                    column2str_length.put("short_name",dataSetHashMap.get(name).getShort_name().length());
                }

                if(column2str_length.containsKey("n_spectra")){
                    if(dataset2n_spectra.get(name) > column2str_length.get("n_spectra")){
                        column2str_length.put("n_spectra",String.valueOf(dataset2n_spectra.get(name)).length());
                    }
                }else{
                    column2str_length.put("n_spectra",String.valueOf(dataset2n_spectra.get(name)).length());
                }

                if(column2str_length.containsKey("n_ms_files")){
                    if(dataset2n_ms_file.get(name) > column2str_length.get("n_ms_files")){
                        column2str_length.put("n_ms_file",String.valueOf(dataset2n_ms_file.get(name)).length());
                    }
                }else{
                    column2str_length.put("n_ms_file",String.valueOf(dataset2n_ms_file.get(name)).length());
                }
            }
        }

        if(column2str_length.containsKey("species") && "species".length() > column2str_length.get("species")){
            column2str_length.put("species","species".length());
        }
        if(column2str_length.containsKey("data_type") && "data_type".length() > column2str_length.get("data_type")){
            column2str_length.put("data_type","data_type".length());
        }
        if(column2str_length.containsKey("short_name") && "short_name".length() > column2str_length.get("short_name")){
            column2str_length.put("short_name","short_name".length());
        }


        if(show_full_information){
            if("n_spectra".length() > column2str_length.get("n_spectra")){
                column2str_length.put("n_spectra","n_spectra".length());
            }
            if("n_ms_file".length() > column2str_length.get("n_ms_file")){
                column2str_length.put("n_ms_file","n_ms_file".length());
            }
            System.out.printf("%s %s %s %s %s %s %s %s %s\n", StringUtils.rightPad("NO.", 3),
                    StringUtils.rightPad("dataset_name", column2str_length.get("dataset_name")),
                    StringUtils.rightPad("short_name", column2str_length.get("short_name")),
                    StringUtils.rightPad("parameter_set", column2str_length.get("parameter_set")),
                    StringUtils.rightPad("species", column2str_length.get("species")),
                    StringUtils.rightPad("data_type", column2str_length.get("data_type")),
                    StringUtils.rightPad("n_spectra", column2str_length.get("n_spectra")),
                    StringUtils.rightPad("n_ms_file", column2str_length.get("n_ms_file")),
                    StringUtils.rightPad("data_link", column2str_length.get("data_link")));


        }else {
            System.out.printf("%s %s %s %s %s %s\n", StringUtils.rightPad("NO.", 3),
                    StringUtils.rightPad("dataset_name", column2str_length.get("dataset_name")),
                    StringUtils.rightPad("parameter_set", column2str_length.get("parameter_set")),
                    StringUtils.rightPad("species", column2str_length.get("species")),
                    StringUtils.rightPad("data_type", column2str_length.get("data_type")),
                    StringUtils.rightPad("data_link", column2str_length.get("data_link")));
        }

        for(String name: names) {
            i = i + 1;
            if (show_full_information) {
                System.out.printf("%s %s %s %s %s %s %s %s %s\n", StringUtils.rightPad(String.valueOf(i), 3),
                        StringUtils.rightPad(name, column2str_length.get("dataset_name")),
                        StringUtils.rightPad(dataSetHashMap.get(name).getShort_name(), column2str_length.get("short_name")),
                        StringUtils.rightPad(dataSetHashMap.get(name).getParameter(), column2str_length.get("parameter_set")),
                        StringUtils.rightPad(dataSetHashMap.get(name).getSpecies(), column2str_length.get("species")),
                        StringUtils.rightPad(dataSetHashMap.get(name).getData_type(), column2str_length.get("data_type")),
                        StringUtils.rightPad(String.valueOf(dataset2n_spectra.get(name)), column2str_length.get("n_spectra")),
                        StringUtils.rightPad(String.valueOf(dataset2n_ms_file.get(name)), column2str_length.get("n_ms_file")),
                        StringUtils.rightPad(dataSetHashMap.get(name).getData_link(), column2str_length.get("data_link")));

            }else{
                System.out.printf("%s %s %s %s %s %s\n", StringUtils.rightPad(String.valueOf(i), 3),
                        StringUtils.rightPad(name, column2str_length.get("dataset_name")),
                        StringUtils.rightPad(dataSetHashMap.get(name).getParameter(), column2str_length.get("parameter_set")),
                        StringUtils.rightPad(dataSetHashMap.get(name).getSpecies(), column2str_length.get("species")),
                        StringUtils.rightPad(dataSetHashMap.get(name).getData_type(), column2str_length.get("data_type")),
                        StringUtils.rightPad(dataSetHashMap.get(name).getData_link(), column2str_length.get("data_link")));
            }
        }

        if(show_full_information){
            long n_spectra = 0L;
            for(String name: dataset2n_spectra.keySet()){
                n_spectra = n_spectra + dataset2n_spectra.get(name);
            }
            long n_ms_file = 0L;
            for(String name: dataset2n_ms_file.keySet()){
                n_ms_file = n_ms_file + dataset2n_ms_file.get(name);
            }
            System.out.println("Total number of datasets:"+dataSetHashMap.size());
            System.out.println("Total number of MS/MS spectra:"+n_spectra);
            System.out.println("Total number of MS file:"+n_ms_file);
        }
    }

    public static HashMap<String, MSDataSet> load_datasets(String dataset_names){
        HashMap<String,MSDataSet> msDataSetHashMap = new HashMap<>();

        if(dataset_names.equalsIgnoreCase("all")){
            msDataSetHashMap = load_default_datasets();
        }else{
            HashMap<String, MSDataSet> default_datasets = load_default_datasets();
            if(dataset_names.contains(",")) {
                // If there is "," in dataset_names, it must be datasets from PepQuery public MS/MS datasets collection.
                String[] names = dataset_names.split(",");
                for (String name : names) {
                    if (default_datasets.containsKey(name)) {
                        msDataSetHashMap.put(name, default_datasets.get(name));
                    } else {
                        Cloger.getInstance().logger.warn("MS/MS dataset is not present:" + name);
                        System.exit(1);
                    }
                }
            }else{
                // could be a single dataset name or a pattern for selecting a set of datasets
                // a single dataset name
                if(default_datasets.containsKey(dataset_names)){
                    msDataSetHashMap.put(dataset_names, default_datasets.get(dataset_names));
                }else{
                    // a pattern for selecting a set of datasets
                    if(dataset_names.length()==1){
                        // selecting data sets based on data_type
                        if(dataset_names.equalsIgnoreCase("w") ||
                                dataset_names.equalsIgnoreCase("g") ||
                                dataset_names.equalsIgnoreCase("a") ||
                                dataset_names.equalsIgnoreCase("p") ||
                                dataset_names.equalsIgnoreCase("u")){

                            for (String name : default_datasets.keySet()) {
                                if (default_datasets.get(name).getData_type().equalsIgnoreCase(dataset_names)) {
                                    msDataSetHashMap.put(name, default_datasets.get(name));
                                }
                            }

                        }else{
                            Cloger.getInstance().logger.error("Invalid value for dataset selection:"+dataset_names);
                            System.exit(1);
                        }


                    }else {
                        Pattern pattern = Pattern.compile(dataset_names);
                        for (String name : default_datasets.keySet()) {
                            Matcher matcher = pattern.matcher(name);
                            if (matcher.find()) {
                                msDataSetHashMap.put(name, default_datasets.get(name));
                            }
                        }
                    }
                }
                if(msDataSetHashMap.size()<=0){
                    Cloger.getInstance().logger.warn("There is no dataset with name or pattern '"+dataset_names+"' in PepQuery public MS/MS datasets collection.");
                    //System.exit(1);
                }
            }
        }

        return msDataSetHashMap;
    }

    public void print(){
        Field[] fields = MSDataSet.class.getDeclaredFields();
        for(Field field: fields){
            try {
                System.out.println(field.getName()+"\t"+field.get(this));
            } catch (IllegalAccessException e) {
                e.printStackTrace();
            }
        }
    }


    /**
     * Get the number of spectra for a dataset
     * @param dataset_name full dataset name
     * @return The number of spectra of a dataset
     * @throws FileNotFoundException
     */
    public long get_the_number_of_spectra(String dataset_name){
        HashMap<String, MSDataSet> dataset = load_datasets(dataset_name);
        MSDataSet msDataSet = dataset.get(dataset_name);
        String ms_summary_file = msDataSet.ms_file + "summary.txt";
        S3Interface.getInstance();
        // total_msms_spectra
        long n_spectra = get_the_number_of_spectra(msDataSet);
        return n_spectra;
    }


    public static long get_the_number_of_spectra(MSDataSet msDataSet){
        String ms_summary_file = msDataSet.ms_file + "summary.txt";
        S3Interface.getInstance();
        // total_msms_spectra
        long n_spectra = 0;
        if (ms_summary_file.startsWith("s3:")) {
            String objPath = ms_summary_file;
            File tempFile = null;
            try {
                tempFile = File.createTempFile(Thread.currentThread().getName(), "_summary.txt");
            } catch (IOException e) {
                e.printStackTrace();
            }
            String tmp_file = tempFile.getAbsolutePath();

            if (S3Interface.getInstance().download(objPath, tmp_file)) {
                try {
                    BufferedReader bReader = new BufferedReader(new FileReader(new File(tmp_file)));
                    String line;
                    while((line = bReader.readLine())!=null){
                        line = line.trim();
                        if(line.startsWith("total_msms_spectra=")){
                            String d[]= line.split("=");
                            n_spectra = Long.parseLong(d[1]);
                        }
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
                tempFile.delete();
            } else {
                Cloger.getInstance().logger.error(objPath + " doesn't exist");
                //System.exit(1);
            }
        }
        return n_spectra;
    }

    public static long get_the_number_of_ms_files(MSDataSet msDataSet){
        String ms_summary_file = msDataSet.ms_file + "summary.txt";
        S3Interface.getInstance();
        long n_ms_files = 0;
        if (ms_summary_file.startsWith("s3:")) {
            String objPath = ms_summary_file;
            File tempFile = null;
            try {
                tempFile = File.createTempFile(Thread.currentThread().getName(), "_summary.txt");
            } catch (IOException e) {
                e.printStackTrace();
            }
            String tmp_file = tempFile.getAbsolutePath();

            if (S3Interface.getInstance().download(objPath, tmp_file)) {
                try {
                    BufferedReader bReader = new BufferedReader(new FileReader(new File(tmp_file)));
                    String line;
                    while((line = bReader.readLine())!=null){
                        line = line.trim();
                        if(line.startsWith("total_ms_file=")){
                            String d[]= line.split("=");
                            n_ms_files = Long.parseLong(d[1]);
                        }
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
                tempFile.delete();
            } else {
                Cloger.getInstance().logger.error(objPath + " doesn't exist");
                //System.exit(1);
            }
        }
        return n_ms_files;
    }
}
