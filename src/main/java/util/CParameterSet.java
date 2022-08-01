package main.java.util;

import com.alibaba.fastjson.JSONArray;
import com.alibaba.fastjson.JSONObject;
import main.java.pg.CParameter;
import org.apache.commons.io.IOUtils;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

public class CParameterSet {
    private int enzyme_index = 1;
    private int max_missed_cleavage = 1;
    private String fixed_modification_index = "1,2,3";
    private String variable_modification_index = "1,2,3";
    private double tol = 10;
    private String tolu = "ppm";
    private double itol = 0.05;
    private String itolu = "Da";

    public void print(){
        String tol_unit;
        if(tolu.equalsIgnoreCase("ppm")){
            tol_unit = "ppm";
        }else{
            tol_unit = "Da";
        }

        String itol_unit;
        if(itolu.equalsIgnoreCase("ppm")){
            itol_unit = "ppm";
        }else{
            itol_unit = "Da";
        }

        StringBuilder sBuilder = new StringBuilder();
        sBuilder.append("\tFixed modification: "+fixed_modification_index+" = "+ CParameter.getModificationString(fixed_modification_index)).append("\n");
        sBuilder.append("\tVariable modification: "+variable_modification_index+" = "+CParameter.getModificationString(variable_modification_index)).append("\n");
        sBuilder.append("\tEnzyme: "+enzyme_index).append("\n");
        sBuilder.append("\tMax Missed cleavages: "+max_missed_cleavage).append("\n");
        sBuilder.append("\tPrecursor mass tolerance: "+tol).append("\n");
        sBuilder.append("\tPrecursor ion mass tolerance unit: "+tol_unit).append("\n");
        sBuilder.append("\tFragment ion mass tolerance: "+itol).append("\n");
        sBuilder.append("\tFragment ion mass tolerance unit: "+itol_unit).append("\n");
        System.out.print(sBuilder.toString());
    }

    public static void print_all(){
        HashMap<String,CParameterSet> parameterSetHashMap = CParameterSet.load_default_parameter_sets();
        ArrayList<String> names = new ArrayList<>(parameterSetHashMap.keySet());
        Collections.sort(names);
        for (String name : names) {
            System.out.println(name);
            parameterSetHashMap.get(name).print();
        }
    }

    public int getEnzyme_index() {
        return enzyme_index;
    }

    public void setEnzyme_index(int enzyme_index) {
        this.enzyme_index = enzyme_index;
    }

    public int getMax_missed_cleavage() {
        return max_missed_cleavage;
    }

    public void setMax_missed_cleavage(int max_missed_cleavage) {
        this.max_missed_cleavage = max_missed_cleavage;
    }

    public String getFixed_modification_index() {
        return fixed_modification_index;
    }

    public void setFixed_modification_index(String fixed_modification_index) {
        this.fixed_modification_index = fixed_modification_index;
    }

    public String getVariable_modification_index() {
        return variable_modification_index;
    }

    public void setVariable_modification_index(String variable_modification_index) {
        this.variable_modification_index = variable_modification_index;
    }

    public double getTol() {
        return tol;
    }

    public void setTol(double tol) {
        this.tol = tol;
    }

    public String getTolu() {
        return tolu;
    }

    public void setTolu(String tolu) {
        this.tolu = tolu;
    }

    public double getItol() {
        return itol;
    }

    public void setItol(double itol) {
        this.itol = itol;
    }

    public String getItolu() {
        return itolu;
    }

    public void setItolu(String itolu) {
        this.itolu = itolu;
    }

    public CParameterSet(){

    }

    public static HashMap<String, CParameterSet> load_default_parameter_sets(){
        InputStream inputStream = MSDataSet.class.getResourceAsStream("/main/resources/search_parameter_set.json");
        String json = null;
        try {
            json = IOUtils.toString(inputStream, StandardCharsets.UTF_8.name());
        } catch (IOException e) {
            e.printStackTrace();
        }
        JSONObject jsonObject = (JSONObject) JSONArray.parse(json);
        HashMap<String,CParameterSet> parameterSetHashMap = new HashMap<>();
        for (Map.Entry<String, Object> entry : jsonObject.entrySet()) {
            //System.out.println(entry.getKey());
            String dataset_id = entry.getKey();
            CParameterSet parameterSet = ((JSONObject) entry.getValue()).toJavaObject(CParameterSet.class);
            parameterSetHashMap.put(dataset_id,parameterSet);
        }
        Cloger.getInstance().logger.info("Load "+parameterSetHashMap.size()+" parameter sets.");
        return parameterSetHashMap;
    }
}
