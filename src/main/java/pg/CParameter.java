package main.java.pg;


import com.compomics.util.experiment.biology.ions.IonFactory;
import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import main.java.PSMMatch.HyperscoreMatch;
import main.java.PSMMatch.JPeptide;
import main.java.db.Database;
import main.java.util.CParameterSet;
import main.java.util.Cloger;
import main.java.util.MSDataSet;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.lang3.StringUtils;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Properties;

/**
 * Parameter class
 */
public class CParameter {

    /**
     * Enum of the types of input target type available.
     */
    public enum TargetInputType {

        //if the target input sequence is not peptide, then use this parameter to specify the type of input event. 1=protein,2=DNA,3=VCF,4=BED,5=GTF;
        peptide("a single peptide or a file contains a list of peptides"),
        protein("a single protein"),
        dna("a dna sequence"),
        vcf("a vcf file"),
        bed("a bed file"),
        gtf("a gtf file");

        /**
         * The description.
         */
        public final String description;

        /**
         * Constructor.
         *
         * @param description the description of the input target type
         */
        TargetInputType(String description) {
            this.description = description;
        }
    }

    /**
     * Enum of the types of search available.
     */
    public enum SearchType {

        novel("novel peptide/protein validation"),
        known("known peptide validation");

        /**
         * The description of the search type.
         */
        public final String description;

        /**
         * Constructor.
         *
         * @param description the description of the search type
         */
        SearchType(String description) {
            this.description = description;
        }
    }

    /**
     * Default p-value cutoff
     */
    public static double p_value_threshold = 0.01;

    /**
     * The cutoff for short peptide
     */
    public static double p_value_threshold_short_peptide = 0.05;

    /**
     * The cutoff for short peptide, length less than this will be treated as short peptide
     */
    public static int short_peptide_length_threshold = 8;

    public static SearchType search_type = SearchType.novel;

    /**
     * Enum of the types of search available.
     */
    public enum ScoreMethod {

        hyper_score("Hyper score"),
        mvh("MVH");

        /**
         * The description of the score method.
         */
        public final String description;

        /**
         * Constructor.
         *
         * @param description the description of the score method
         */
        ScoreMethod(String description) {
            this.description = description;
        }
    }

    /**
     * Scoring algorithm:
     * 1: Hyperscore
     * 2: MVH
     */
    public static int scoreMethod = 1;


    public static TargetInputType input_target_type = TargetInputType.peptide;
    public static String input_target_file = "";
    public static String dna_frame = "0";
    public static String annotation_data_dir = "";

    public static boolean use_tag_file = false;
    public static String tag_file = "";

    public static boolean do_modification_search_validation = true;
    public static boolean plot_spectrum = false;
    public static boolean generate_psm_annotation_data = false;

    public static boolean add_extra_validation = false;

    /**
     * Reference protein database
     */
    public static String db = "";

    /**
     * internal reference protein database
     * This is usually generated in validating known protein.
     * Default value must be "".
     */
    public static String internal_db = "";


    public static String cmd = "-";

    public static double tol = 10;

    /**
     * decoy protein search
     */
    public static boolean add_decoy = false;

    /**
     * reverse or random
     */
    public static String decoy_method = "reverse";

    /**
     * true: ppm, false:da
     */
    public static String tolu = "ppm";

    /**
     * tolerance for MS2, the unit is Da
     */
    public static double itol = 0.6;

    /**
     * true: ppm, false:da
     */
    public static boolean itolu = false;

    public static int minPeptideLength = 7;

    public static int maxPeptideLength = 45;

    public static double minPeptideMass = 500.0;

    public static double maxPeptideMass = 10000.0;

    /**
     * The maximum number of variable modifications
     */
    public static int maxVarMods = 3;

    /**
     * The maximum number of allowed modifications (variable and fixed) on the same amino acid
     * Default only one modification (either variable or fixed) is allowed on the same position
     */
    public static int maxModsPerAA = 1;

    /**
     * The maximum number of allowed missed cleavage sites
     */
    public static int maxMissedCleavages = 2;

    /**
     * options.addOption("e",true,"1:Trypsin (default), 2:Trypsin (no P rule), 3:Arg-C, 4:Arg-C (no P rule), 5:Arg-N, 6:Glu-C, 7:Lys-C");
     */
    public static int enzyme = 1;

    /**
     * The maximum number of random peptides for p-value calculation, default is 10000
     */
    public static int nRandomPeptides = 10000;


    /**
     * The min score
     */
    public static double minScore = 12.0;

    /**
     * The minimum number of matched peaks for scoring for target peptides.
     */
    public static int min_n_matched_peaks = 4;

    /**
     * The min number of peaks for an MS/MS spectrum
     */
    public static int minPeaks = 10;


    /**
     * 1: CID/HCD (default), 2: ETD
     */
    public static int fragmentMethod = 1;


    /**
     * When an MS/MS spectrum doesn't have precursor change information, the min charge state will be considered.
     */
    public static int minCharge = 2;
    /**
     * When an MS/MS spectrum doesn't have precursor change information, the max charge state will be considered.
     */
    public static int maxCharge = 3;


    public static int cpu = 0;


    public static double intensityThreshold = 0.03;


    /**
     * Fixed modification
     */
    public static String fixMods = "1";

    /**
     * Variable modification
     */
    public static String varMods = "2";

    /**
     * This parameter is used to control whether add the AA substitution modifications when performing the modification
     * filtering. In default, it's false. But when performing missing protein identification, we can set it as true.
     */
    public static boolean addAAsubstitutionMods = false;

    /**
     * When perform unrestricted modification searching, how to filter the result.
     * When the value is true, then filtering is score(ptm) >= score(target peptide).
     * When the value is false, then filtering is score(ptm) > score(target peptide). This is the default value.
     * We recommend to set it as true when require very stringent filtering, for example performing missing protein
     * identification.
     */
    public static boolean unrestrictedSearchWithEqualAndBetterScore = false;

    /**
     * The number of PSMs with better score needed to be reported in unrestricted modification searching-based
     * filtering.
     */
    public static boolean fast_model = false;


    public final static int PepMappingType_FM = 0;

    /**
     * For example, 1,2,3
     */
    public static String isotope_error = "0";
    private final static double isotope_error_mass = 1.00335;

    public static ArrayList<Double> get_peptide_mass_with_isotope_error(double pep_mass){
        ArrayList<Double> mass = new ArrayList<>();
        // It must contain 0;
        // For example isotope_error = "-1,0,1,2"
        // experiment precursor mass = 50, then theoretical mass should be 51,50,49,48
        // if peptide mass is 50, then -1 error, precursor mass is 49
        // if peptide mass is 50, then 0 error, precursor mass is 50
        // if peptide mass is 50, then 1 error, precursor mass is 51
        // if peptide mass is 50, then 2 error, precursor mass is 52
        if(isotope_error.contains("0")) {
            String[] ie = isotope_error.split(",");
            for(String e:ie){
                double new_mass = pep_mass + Double.parseDouble(e)*isotope_error_mass;
                mass.add(new_mass);
                //System.out.println(pep_mass+"\t"+ new_mass);
            }
        }else{
            Cloger.getInstance().logger.error("Invalid isotope error input:"+isotope_error);
            System.exit(1);
        }
        return mass;
    }

    /**
     * Calculate mass tolerance
     * @param spectrum A Spectrum object
     * @param peptide A Peptide object
     * @return An array: tol in ppm, tol in Da and isotope error.
     */
    public static ArrayList<Double> get_mass_error(Spectrum spectrum, Peptide peptide){
        double exp_mass = spectrum.getPrecursor().getMass(spectrum.getPrecursor().possibleCharges[0]);
        double pep_mass = JPeptide.getMass(peptide);
        ArrayList<Double> tol_res = CParameter.get_mass_error(exp_mass,pep_mass);
        return tol_res;
    }

    /**
     * Calculate the delta mass
     * @param precursor_mass spectrum precursor mass
     * @param pep_mass peptide mass
     * @param tol_u the unit of delta mass, ppm or da
     * @return The delta mass in the unit of "tol_u" (ppm or da)
     */
    public static double get_mass_error(double precursor_mass, double pep_mass, String tol_u){
        double delta_mass;
        if(tol_u.equalsIgnoreCase("ppm")){
            delta_mass = (1.0e6) * (pep_mass - precursor_mass) / pep_mass;
        }else{
            delta_mass = pep_mass - precursor_mass;
        }
        return delta_mass;
    }

    /**
     * Calculate mass tolerance
     * @param precursor_mass The precursor mass for a spectrum
     * @param pep_mass The peptide mass
     * @return An array: tol in ppm, tol in Da and isotope error.
     */
    public static ArrayList<Double> get_mass_error(double precursor_mass, double pep_mass){
        ArrayList<Double> res = new ArrayList<>(3);
        boolean tol_valid = false;
        if(tolu.equalsIgnoreCase("ppm")) {
            if(isotope_error.equalsIgnoreCase("0")) {
                double tolPpm = (1.0e6) * (pep_mass - precursor_mass) / pep_mass;
                if(Math.abs(tolPpm) <= CParameter.tol){
                    double tolDa = pep_mass - precursor_mass;
                    res.add(tolPpm);
                    res.add(tolDa);
                    res.add(0.0);
                    tol_valid = true;
                }
            }else{
                String[] ies = isotope_error.split(",");
                for(String ie:ies){
                    double i_ie = Integer.parseInt(ie);
                    // for example, if isotope error is 1 and precursor mass is 100, then the mono mass of precursor
                    // is 100-1=99
                    double tolPpm = (1.0e6) * (pep_mass - (precursor_mass - i_ie*isotope_error_mass)) / pep_mass;
                    if(Math.abs(tolPpm) <= CParameter.tol){
                        double tolDa = pep_mass - (precursor_mass - i_ie*isotope_error_mass);
                        res.add(tolPpm);
                        res.add(tolDa);
                        res.add(i_ie);
                        tol_valid = true;
                        break;
                    }
                }
            }
        }else{
            // da
            if(isotope_error.equalsIgnoreCase("0")) {
                double tolDa = pep_mass - precursor_mass;
                double tolPpm = (1.0e6) * (pep_mass - precursor_mass) / pep_mass;
                if(Math.abs(tolDa) <= CParameter.tol){
                    res.add(tolPpm);
                    res.add(tolDa);
                    res.add(0.0);
                    tol_valid = true;
                }
            }else{
                String[] ies = isotope_error.split(",");
                for(String ie:ies){
                    double i_ie = Integer.parseInt(ie);
                    double tolPpm = (1.0e6) * (pep_mass - (precursor_mass - i_ie*isotope_error_mass)) / pep_mass;
                    double tolDa = pep_mass - (precursor_mass - i_ie*isotope_error_mass);
                    if(Math.abs(tolDa) <= CParameter.tol){
                        res.add(tolPpm);
                        res.add(tolDa);
                        res.add(i_ie);
                        tol_valid = true;
                        break;
                    }
                }
            }
        }
        /**
        if(!tol_valid){
            Cloger.getInstance().logger.error("Mass tol calculation error: precursor_mass="+precursor_mass+", pep_mass="+pep_mass);
            //System.exit(1);
        }**/
        return res;
    }


    public static ArrayList<Double> get_precursor_mass_with_isotope_error(double precursor_mass){
        ArrayList<Double> mass = new ArrayList<>();
        // It must contain 0;
        // For example isotope_error = "-1,0,1,2"
        // experiment precursor mass = 50, then theoretical mass should be 51,50,49,48
        // if peptide mass is 50, then -1 error, precursor mass is 49
        // if peptide mass is 50, then 0 error, precursor mass is 50
        // if peptide mass is 50, then 1 error, precursor mass is 51
        // if peptide mass is 50, then 2 error, precursor mass is 52
        if(isotope_error.contains("0")) {
            if(isotope_error.equals("0")) {
                mass.add(precursor_mass);
            }else{
                String[] ie = isotope_error.split(",");
                for (String e : ie) {
                    double new_mass = precursor_mass - Double.parseDouble(e) * isotope_error_mass;
                    mass.add(new_mass);
                    //System.out.println(precursor_mass+"\t"+ new_mass);
                }
            }
        }else{
            Cloger.getInstance().logger.error("Invalid isotope error input:"+isotope_error);
            System.exit(1);
        }
        return mass;
    }

    /**
     * Peptide mapping type. Use FM index method.
     */
    public static int peptideMappingType = PepMappingType_FM;

    /**
     * Output folder
     */
    public static String outdir = "./";

    /**
     * 1: index (1-based), 2: spectrum title in MGF file
     */
    public static int indexType = 1;


    public static void print() throws IOException {

        String tol_unit;
        if(tolu.equalsIgnoreCase("ppm")){
            tol_unit = "ppm";
        }else{
            tol_unit = "Da";
        }

        String itol_unit;
        if(itolu){
            itol_unit = "ppm";
        }else{
            itol_unit = "Da";
        }

        String score_name;
        if(scoreMethod == 1){
            score_name = "Hyperscore";
        }else if(scoreMethod == 2){
            score_name = "MVH";
        }else{
            score_name = "Invalid";
        }

        System.out.println("#############################################");
        System.out.println("PepQuery parameter:");

        StringBuilder sBuilder = new StringBuilder();
        sBuilder.append("PepQuery version: " + getVersion()).append("\n");
        if(!cmd.equalsIgnoreCase("-")){
            sBuilder.append("PepQuery command line: " + cmd).append("\n");
        }
        sBuilder.append("Fixed modification: "+fixMods+" = "+getModificationString(fixMods)).append("\n");
        sBuilder.append("Variable modification: "+varMods+" = "+getModificationString(varMods)).append("\n");
        sBuilder.append("Max allowed variable modification: "+maxVarMods).append("\n");
        sBuilder.append("Add AA substitution: "+addAAsubstitutionMods).append("\n");
        sBuilder.append("Enzyme: "+enzyme+" = "+DatabaseInput.getEnzymeByIndex(CParameter.enzyme).getName()).append("\n");
        sBuilder.append("Max Missed cleavages: "+maxMissedCleavages).append("\n");
        sBuilder.append("Precursor mass tolerance: "+tol).append("\n");
        sBuilder.append("Range of allowed isotope peak errors: "+isotope_error).append("\n");


        sBuilder.append("Precursor ion mass tolerance unit: "+tol_unit).append("\n");
        sBuilder.append("Fragment ion mass tolerance: "+itol).append("\n");
        sBuilder.append("Fragment ion mass tolerance unit: "+itol_unit).append("\n");
        sBuilder.append("Scoring algorithm: "+scoreMethod+" = "+score_name).append("\n");
        sBuilder.append("Min score: "+minScore).append("\n");
        sBuilder.append("Min peaks: "+minPeaks).append("\n");

        sBuilder.append("Min peptide length: "+minPeptideLength).append("\n");
        sBuilder.append("Max peptide length: "+maxPeptideLength).append("\n");

        sBuilder.append("Min peptide mass: "+minPeptideMass).append("\n");
        sBuilder.append("Max peptide mass: "+maxPeptideMass).append("\n");

        sBuilder.append("Random peptide number: "+nRandomPeptides).append("\n");
        sBuilder.append("Fast mode: "+fast_model).append("\n");
        sBuilder.append("CPU: "+cpu).append("\n");
        System.out.print(sBuilder.toString());
        System.out.println("#############################################");


        // Save parameters and command line information to a file
        String para_file = outdir + "/parameter.txt";
        BufferedWriter bWriter = new BufferedWriter(new FileWriter(new File(para_file)));
        //bWriter.write("parameter\tvalue\n");
        bWriter.write(sBuilder.toString());
        bWriter.close();



    }


    public static void init(){
        HyperscoreMatch.generateFactorialValues(60);
        // When multiple threads access this, it may cause some problems.
        /**
         * Exception in thread "pool-4-thread-19" java.util.ConcurrentModificationException
         * 	at java.base/java.util.HashMap$HashIterator.nextNode(HashMap.java:1493)
         * 	at java.base/java.util.HashMap$KeyIterator.next(HashMap.java:1516)
         * 	at java.base/java.util.AbstractCollection.addAll(AbstractCollection.java:351)
         * 	at java.base/java.util.HashSet.<init>(HashSet.java:120)
         * 	at com.compomics.util.experiment.biology.ions.IonFactory.getFragmentIons(IonFactory.java:231)
         *
         *    public static HashSet<String> getDefaultNeutralLosses() {
         *         if (defaultNeutralLosses == null) {
         *             setDefaultNeutralLosses();
         *         }
         *         return defaultNeutralLosses;
         *     }
         *
         *         private static synchronized void setDefaultNeutralLosses() {
         *         defaultNeutralLosses = new HashSet<>(2);
         *         defaultNeutralLosses.add(NeutralLoss.H2O.name);
         *         defaultNeutralLosses.add(NeutralLoss.NH3.name);
         *     }
         *
         */
        IonFactory.getDefaultNeutralLosses();
    }

    public static String getModificationString (String modValue){
        ModificationGear.getInstance();
        ArrayList<String> mods = new ArrayList<>();
        if(!modValue.equalsIgnoreCase("0") && modValue != null && !modValue.isEmpty()){
            String fm[] = modValue.split(",");
            for(int i=0;i<fm.length;i++){
                String ptMname = ModificationGear.getInstance().getPTMname(Integer.valueOf(fm[i]));
                mods.add(ptMname);
            }
        }
        String res = "-";
        if(mods.size()>=1){
            res = StringUtils.join(mods,',');
        }
        return(res);

    }

    /**
     * Get version
     * @return
     */
    public static String getVersion(){
        Properties properties = new Properties();
        try {
            properties.load(CParameter.class.getClassLoader().getResourceAsStream("project.properties"));
        } catch (IOException e) {
            e.printStackTrace();
        }
        return properties.getProperty("version");
    }

    /**
     * Update parameters in CParameter using msmsDataSet
     * @param msDataSet MSDataSet
     */
    public static void updateCParameter(MSDataSet msDataSet){
        CParameter.enzyme = msDataSet.parameterSet.getEnzyme_index();
        CParameter.maxMissedCleavages = msDataSet.parameterSet.getMax_missed_cleavage();
        CParameter.do_modification_search_validation = true;
        CParameter.itol = msDataSet.parameterSet.getItol();
        CParameter.tol = msDataSet.parameterSet.getTol();
        CParameter.tolu = msDataSet.parameterSet.getTolu();
        CParameter.fixMods = msDataSet.parameterSet.getFixed_modification_index();
        CParameter.varMods = msDataSet.parameterSet.getVariable_modification_index();

    }


    public static void updateCParameter(CParameterSet parameterSet){
        CParameter.enzyme = parameterSet.getEnzyme_index();
        CParameter.maxMissedCleavages = parameterSet.getMax_missed_cleavage();
        CParameter.itol = parameterSet.getItol();
        CParameter.tol = parameterSet.getTol();
        CParameter.tolu = parameterSet.getTolu();
        CParameter.fixMods = parameterSet.getFixed_modification_index();
        CParameter.varMods = parameterSet.getVariable_modification_index();
    }

    /**
     * Update parameters in CParameter using command line options
     * @param cmd CommandLine object
     */
    public static void updateCParameter(CommandLine cmd){

        if(cmd.hasOption("p")){
            String parameter_set_name = cmd.getOptionValue("p");
            HashMap<String,CParameterSet> parameterSetHashMap = CParameterSet.load_default_parameter_sets();
            if(parameterSetHashMap.containsKey(parameter_set_name)){
                Cloger.getInstance().logger.info("Update parameter using parameter set:"+parameter_set_name);
                //parameterSetHashMap.get(parameter_set_name).print();
                CParameter.updateCParameter(parameterSetHashMap.get(parameter_set_name));
            }else{
                Cloger.getInstance().logger.error("Invalid parameter name:"+parameter_set_name);
                System.exit(1);
            }


        }
        // don't set default value in this function.
        // only update a parameter if it's provided a command line option.
        if(cmd.hasOption("db")){
            if(CParameter.db.isEmpty()) {
                if(cmd.getOptionValue("db").contains(":")){
                    String db_p[] = cmd.getOptionValue("db").split(":");
                    String odir = "." + File.separator + "database";
                    if(cmd.hasOption("o")){
                        odir = cmd.getOptionValue("o") + File.separator + "database";
                    }
                    String db = Database.download_db(db_p[0],db_p[1],"",odir);
                    CParameter.db = db;
                }else {
                    CParameter.db = cmd.getOptionValue("db");
                }
                CParameter.internal_db = CParameter.db;
            }
        }else{
            Cloger.getInstance().logger.error("Protein reference database is not provided!");
            System.exit(1);
        }

        if(cmd.hasOption("fast")){
            CParameter.fast_model = true;
        }

        if(cmd.hasOption("e")){
            CParameter.enzyme = Integer.parseInt(cmd.getOptionValue("e"));
        }
        if(cmd.hasOption("c")){
            int missCleavages = Integer.parseInt(cmd.getOptionValue("c"));
            if(missCleavages<0){
                missCleavages = 0;
            }
            CParameter.maxMissedCleavages = missCleavages;
        }

        if(cmd.hasOption("ti")){
            CParameter.isotope_error = cmd.getOptionValue("ti");
        }

        if(CParameter.enzyme == 0){
            CParameter.maxMissedCleavages = 100;
        }

        if(cmd.hasOption("tol")){
            CParameter.tol = Double.parseDouble(cmd.getOptionValue("tol"));
        }

        if(cmd.hasOption("tolu")){
            if(cmd.getOptionValue("tolu").equalsIgnoreCase("ppm")){
                CParameter.tolu = "ppm";
            }else{
                CParameter.tolu = "da";
            }
        }


        if(cmd.hasOption("itol")){
            CParameter.itol = Double.parseDouble(cmd.getOptionValue("itol"));
        }

        //JMatch.setMs2tol(itol);

        //int method = 1; // hyperscore
        if(cmd.hasOption("m")){
            CParameter.scoreMethod = Integer.parseInt(cmd.getOptionValue("m"));
        }

        //int nRandomPeptides = 1000;
        if(cmd.hasOption("n")){
            CParameter.nRandomPeptides = Integer.parseInt(cmd.getOptionValue("n"));
        }


        //int cpu = 1;
        if(cmd.hasOption("cpu")){
            CParameter.cpu = Integer.parseInt(cmd.getOptionValue("cpu"));
            if(CParameter.cpu==0){
                CParameter.cpu = Runtime.getRuntime().availableProcessors();
            }
        }else{
            CParameter.cpu = Runtime.getRuntime().availableProcessors();
        }


        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if(cmd.hasOption("minScore")){
            CParameter.minScore = Double.parseDouble(cmd.getOptionValue("minScore"));
        }

        if(cmd.hasOption("minPeaks")){
            CParameter.minPeaks = Integer.parseInt(cmd.getOptionValue("minPeaks"));
        }

        //double limit = 0.03;
        if(cmd.hasOption("limit")){
            CParameter.intensityThreshold = Double.parseDouble(cmd.getOptionValue("limit"));
        }

        // int maxVar = 3;
        if(cmd.hasOption("maxVar")){
            CParameter.maxVarMods = Integer.parseInt(cmd.getOptionValue("maxVar"));
        }

        //int minPeptideLength = 7;
        if(cmd.hasOption("minLength")){
            CParameter.minPeptideLength = Integer.parseInt(cmd.getOptionValue("minLength"));
        }

        //int maxPeptideLength = 45;
        if(cmd.hasOption("maxLength")){
            CParameter.maxPeptideLength = Integer.parseInt(cmd.getOptionValue("maxLength"));
        }

        //int minCharge = 2;
        if(cmd.hasOption("minCharge")){
            CParameter.minCharge = Integer.parseInt(cmd.getOptionValue("minCharge"));
        }

        //int maxCharge = 3;
        if(cmd.hasOption("maxCharge")){
            CParameter.maxCharge = Integer.parseInt(cmd.getOptionValue("maxCharge"));
        }

        //String fixModOptionValue = "6";
        //String varModOptionValue = "107";

        if(cmd.hasOption("fixMod")){
            CParameter.fixMods = cmd.getOptionValue("fixMod");
            // no fixed modification
            if(CParameter.fixMods.equalsIgnoreCase("") ||
                    CParameter.fixMods.equalsIgnoreCase("no") ||
                    CParameter.fixMods.equalsIgnoreCase("0")){
                CParameter.fixMods = "0";
            }
        }
        if(cmd.hasOption("varMod")){
            CParameter.varMods = cmd.getOptionValue("varMod");
            if(CParameter.varMods.equalsIgnoreCase("") ||
                    CParameter.varMods.equalsIgnoreCase("no") ||
                    CParameter.varMods.equalsIgnoreCase("0")){
                CParameter.varMods = "0";
            }
        }

        // For most of the applications, we don't need to consider these modifications in
        // unrestricted modification-based filtering.
        if(cmd.hasOption("aa")){
            CParameter.addAAsubstitutionMods = true;
        }

        if(cmd.hasOption("hc")){
            CParameter.unrestrictedSearchWithEqualAndBetterScore = true;
        }

        if(cmd.hasOption("indexType")){
            CParameter.indexType = Integer.parseInt(cmd.getOptionValue("indexType"));
        }

        if(cmd.hasOption("i")){
            CParameter.input_target_file = cmd.getOptionValue("i");
        }else{
            Cloger.getInstance().logger.error("Input for -i is required!");
            System.exit(1);
        }

        if(cmd.hasOption("s")){
            // Perform target protein identification. In this case, the value for "-i" is a protein ID.
            if(cmd.getOptionValue("s").equalsIgnoreCase("1") || cmd.getOptionValue("s").equalsIgnoreCase("novel")){
                // novel peptide/protein validation
                CParameter.search_type = SearchType.novel;
            }else if(cmd.getOptionValue("s").equalsIgnoreCase("2") || cmd.getOptionValue("s").equalsIgnoreCase("known")){
                CParameter.search_type = CParameter.SearchType.known;
            }else{
                Cloger.getInstance().logger.warn("The value for -'s' is not valid:"+cmd.getOptionValue("s"));
                System.exit(1);
            }
        }else{
            CParameter.search_type = SearchType.novel;
        }

        if(cmd.hasOption("decoy")){
            CParameter.add_decoy = true;
        }

        if(cmd.hasOption("t")){
            if(cmd.getOptionValue("t").equalsIgnoreCase("peptide") || cmd.getOptionValue("t").equalsIgnoreCase("pep")) {
                CParameter.input_target_type = CParameter.TargetInputType.peptide;
                if(CParameter.search_type.equals(SearchType.novel)){
                    Cloger.getInstance().logger.info("Task type: novel peptide identification");
                }else{
                    Cloger.getInstance().logger.info("Task type: known peptide identification");
                }
            }else if(cmd.getOptionValue("t").equalsIgnoreCase("protein") || cmd.getOptionValue("t").equalsIgnoreCase("pro")){
                CParameter.input_target_type = CParameter.TargetInputType.protein;
                if(CParameter.search_type.equals(SearchType.novel)){
                    Cloger.getInstance().logger.info("Task type: novel protein identification");
                }else{
                    Cloger.getInstance().logger.info("Task type: known protein identification");
                }
            }else if(cmd.getOptionValue("t").equalsIgnoreCase("dna")){
                CParameter.input_target_type = CParameter.TargetInputType.dna;
                CParameter.search_type = SearchType.novel;
                Cloger.getInstance().logger.info("Task type: novel DNA sequence identification");
            }else if(cmd.getOptionValue("t").equalsIgnoreCase("vcf")){
                CParameter.input_target_type = CParameter.TargetInputType.vcf;
                CParameter.search_type = SearchType.novel;
                Cloger.getInstance().logger.info("Task type: variant peptide (vcf) identification");
            }else if(cmd.getOptionValue("t").equalsIgnoreCase("bed")){
                CParameter.input_target_type = CParameter.TargetInputType.bed;
                CParameter.search_type = SearchType.novel;
                Cloger.getInstance().logger.info("Task type: novel peptide (bed) identification");
            }else if(cmd.getOptionValue("t").equalsIgnoreCase("gtf")){
                CParameter.input_target_type = CParameter.TargetInputType.gtf;
                CParameter.search_type = SearchType.novel;
                Cloger.getInstance().logger.info("Task type: novel peptide (gtf) identification");
            }else{
                Cloger.getInstance().logger.error("Invalid input type for -t: "+cmd.getOptionValue("t"));
                System.exit(1);
            }
        }else{
            CParameter.input_target_type = CParameter.TargetInputType.peptide;
            if(CParameter.search_type.equals(SearchType.novel)){
                Cloger.getInstance().logger.info("Task type: novel peptide identification");
            }else{
                Cloger.getInstance().logger.info("Task type: known peptide identification");
            }
        }
        if(cmd.hasOption("frame")){
            CParameter.dna_frame = cmd.getOptionValue("frame");
        }
        if(cmd.hasOption("anno")) {
            CParameter.annotation_data_dir = cmd.getOptionValue("anno");
        }



        //if(cmd.hasOption("um")) {
        //    CParameter.do_modification_search_validation = true;
        //}

        if(cmd.hasOption("plot")){
            //CParameter.plot_spectrum = true;
            CParameter.generate_psm_annotation_data = true;
        }

        if(cmd.hasOption("x")){
            CParameter.add_extra_validation = true;
        }
    }

    public static boolean passed_p_value_validation(int peptide_length, double p_value){
        if(peptide_length <= CParameter.short_peptide_length_threshold &&
                p_value <= CParameter.p_value_threshold_short_peptide){
            return true;
        }else if(peptide_length > CParameter.short_peptide_length_threshold &&
                p_value <= CParameter.p_value_threshold){
            return true;
        }else{
            return false;
        }
    }

    public static boolean passed_psm_validation(int peptide_length, double p_value, int n_db, int rank){

        boolean p_value_passed;
        if(peptide_length <= CParameter.short_peptide_length_threshold &&
                p_value <= CParameter.p_value_threshold_short_peptide){
            p_value_passed =  true;
        }else if(peptide_length > CParameter.short_peptide_length_threshold &&
                p_value <= CParameter.p_value_threshold){
            p_value_passed =  true;
        }else{
            p_value_passed =  false;
        }

        if(p_value_passed && n_db == 0 && rank == 1){
            return true;
        }else{
            return false;
        }
    }

    public static boolean passed_psm_validation(int peptide_length, double p_value, int n_db, int rank, int n_ptm){
        boolean passed = passed_psm_validation(peptide_length, p_value, n_db, rank);
        if(passed && n_ptm ==0){
            return true;
        }else{
            return false;
        }
    }




}
