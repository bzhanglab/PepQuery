package main.java.pg;

import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.io.mass_spectrometry.mgf.MgfFileIterator;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import com.compomics.util.gui.waiting.waitinghandlers.WaitingHandlerCLIImpl;
import com.compomics.util.waiting.WaitingHandler;
import main.java.OpenModificationSearch.ModificationDB;
import main.java.download.FileDownload;
import main.java.index.BuildMSlibrary;
import main.java.index.IndexWorker;
import main.java.plot.GeneratePSMAnnotation;
import main.java.plot.PeakAnnotation;
import main.java.PSMMatch.*;
import main.java.util.CParameterSet;
import main.java.util.Cloger;
import main.java.util.MSDataSet;
import org.apache.commons.cli.*;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import java.io.*;
import java.sql.SQLException;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;



public class PeptideSearchMT {

    public static boolean debug = false;


    public static void main(String[] args) throws ParseException, IOException, InterruptedException, SQLException, ClassNotFoundException {


        Options options = new Options();
        options.addOption("i", true, "The value for this parameter could be (1) a single peptide sequence; " +
                "(2) multiple peptide sequences separated by ','; " +
                "(3) a file contains peptide sequence(s). If it is a file contains peptide sequence(s), each row is a single peptide sequence and there is no header in the file; " +
                "(4) a file contains peptide sequence and spectrum ID and they are separated by tab '\\t'. There is no header in the file. This is used when perform PSM level validation; " +
                "(5) a protein sequence. This is used when perform novel protein identification; " +
                "(6) a protein ID. This is used for targeted known protein identification; " +
                "(7) a DNA sequence; " +
                "(8) a VCF file; "+
                "(9) a bed file; "+
                "(10) a GTF file.");
        options.addOption("t",true,"Input type for parameter -i: peptide (or pep) when the setting for -i belongs to (1)-(4), protein (or pro) when the setting for -i belongs to (5)-(6), DNA (or dna), VCF (or vcf), BED (or bed), GTF or (gtf). The default is peptide.");
        options.addOption("s", true, "Task type: 1 => novel peptide/protein validation (default), 2 => known peptide/protein validation. Default is 1.");
        options.addOption("ms", true, "MS/MS data used for identification. MGF, mzML, mzXML, mgf.gz, mzML.gz, mzXML.gz and raw formats are supported. USI is also supported.");
        // For MGF file, what type of index should be used.
        options.addOption("indexType", true, "When the input MS/MS data is in MGF format, what type of index will be used: 1 => index (1-based), 2 => spectrum title in MGF file. Default is 1.");
        options.addOption("b", true, "MS/MS dataset ID(s) to search. Two types of dataset ID are supported: the dataset ID from indexed PepQueryDB and dataset ID from public proteomics data repositories including PRIDE, MassIVE, jPOSTrepo and iProX. " +
                "Multiple datasets from PepQueryDB must be separated by comma. A pattern to match datasets in PepQueryDB is also supported, for example, use '-b CPTAC' to search all datasets contain 'CPTAC'. " +
                "In addition, dataset selection from PepQueryDB based on data type (w:global proteome, p:phosphorylation, a:acetylation, u:ubiquitination, g:glycosylation) is also supported. For example, use '-b p' to search all phosphoproteomics datasets in PepQueryDB. "+
                "Use '-b show' to present all MS/MS datasets available in PepQueryDB. " +
                "For dataset ID from public proteomics data repositories, one dataset is supported for each analysis. For example, use '-b PXD000529' to use all MS/MS data from dataset PXD000529 or " +
                "use '-b PXD000529:LM3' to use data files containing LM3 from dataset PXD000529.");
        //options.addOption("pep", true, "Peptide sequence which you want to search");
        options.addOption("db", true, "A protein reference database in FASTA format. For example, a protein sequence database from RefSeq, GENCODE, Ensembl or UniProt. A string like swissprot:human, refseq:human or gencode:human is also accepted. If later is the case, PepQuery will automatically download protein database from online databases. ");
        options.addOption("frame",true,"The frame to translate DNA sequence to protein. The right format is like this: '1,2,3,4,5,6', '1,2,3', '1' or '0' which means to keep the longest frame. In default, for each frame only the longest protein is used. It is only used when the input for '-i' is a DNA sequence.");
        options.addOption("anno",true,"Annotation files folder for VCF/BED/GTF");
        options.addOption("decoy",false,"In known protein identification mode, try to identity the decoy version of the selected target protein. Default is false.");

        // MS/MS data generation related parameters
        options.addOption("e",true,"Enzyme used for protein digestion. 0:Non enzyme, 1:Trypsin (default), 2:Trypsin (no P rule), 3:Arg-C, 4:Arg-C (no P rule), 5:Arg-N, 6:Glu-C, 7:Lys-C");
        options.addOption("c",true,"The max missed cleavages, default is 2");
        options.addOption("minLength", true, "The minimum length of peptide to consider, default is 7");
        options.addOption("maxLength", true, "The maximum length of peptide to consider, default is 45");
        options.addOption("tol", true, "Precursor ion m/z tolerance, default is 10");
        options.addOption("tolu", true, "The unit of precursor ion m/z tolerance, default is ppm");
        options.addOption("itol", true, "Fragment ion m/z tolerance in Da, default is 0.6");
        options.addOption("fragmentMethod", true, "1: CID/HCD (default), 2: ETD");
        options.addOption("fixMod",true,"Fixed modification, the format is like : 1,2,3. Use '-printPTM' to show all supported modifications. Default is 1 (Carbamidomethylation(C)[57.02]). " +
                "If there is no fixed modification, set it as '-fixMod no' or '-fixMod 0'.");
        options.addOption("varMod",true,"Variable modification, the format is same with -fixMod. Default is 2 (Oxidation(M)[15.99]). "+
                "If there is no variable modification, set it as '-varMod no' or '-varMod 0'.");
        options.addOption("maxVar",true,"Max number of variable modifications, default is 3");
        options.addOption("aa",false,"Whether or not to consider aa substitution modifications when perform modification filtering. In default, don't consider.");
        options.addOption("printPTM",false,"Print PTMs");
        //options.addOption("mv",true,"Max number of variable modifications, default is 3");
        options.addOption("ti",true,"Range of allowed isotope peak errors, such as '0,1'; Default: 0");
        options.addOption("p",true,"MS/MS searching parameter set name. Use '-p show' to show all predefined parameter sets.");


        // scoring related parameters
        options.addOption("m",true,"Scoring method: 1=HyperScore (default), 2=MVH");
        // When perform unrestricted modification searching, how to filter the result.
        // When the value is true, then filtering is score(ptm) >= score(target peptide).
        // When the value is false, then filtering is score(ptm) > score(target peptide). This is the default value.
        // We recommend to set it as true when require very stringent filtering,
        // for example performing missing protein identification.
        options.addOption("hc",false,"When perform validation with unrestricted modification searching (UMS), whether or not to use more stringent criterion. TRUE: score(UMS)>=score(targetPSM); FALSE: score(UMS)>score(targetPSM), default");
        options.addOption("x",false,"Add extra score validation. It means use two scoring algorithms for peptide identification.");
        options.addOption("fast",false,"Choose to use the fast mode for searching or not. In fast mode, only one better match from reference peptide-based competitive filtering steps will be returned. A peptide identified or not is not affected by this setting.");
        options.addOption("n",true,"The number of random peptides generated for p-value calculation, default is 10000.");
        // spectra pre-processing
        options.addOption("minPeaks",true,"Min peaks in spectrum, default is 10");
        // In general, the larger the better for the score
        options.addOption("minScore",true,"Minimum score to consider for peptide searching, default is 12");

        // output
        options.addOption("o", true, "Output directory");

        // Unrestricted modification peptide identification
        //options.addOption("um",false,"Validation with unrestricted modification searching");

        options.addOption("cpu",true,"The number of cpus used, default is to use all available CPUs");

        options.addOption("minCharge", true, "The minimum charge to consider if the charge state is not available, default is 2");
        options.addOption("maxCharge", true, "The maximum charge to consider if the charge state is not available, default is 3");

        options.addOption("dc",true,"Tool path for Raw MS/MS data conversion: for example, /home/test/ThermoRawFileParser/ThermoRawFileParser.exe");
        options.addOption("plot", false, "Generate PSM annotation data for plot.");

        options.addOption("h", false, "Help");

        CommandLineParser parser = new DefaultParser(false);
        CommandLine cmd = parser.parse(options, args);
        if (cmd.hasOption("h") || cmd.hasOption("help") || args.length == 0) {
            HelpFormatter f = new HelpFormatter();
            f.setWidth(100);
            f.setOptionComparator(null);
            System.out.println("java -Xmx2G -jar pepquery.jar");
            f.printHelp("Options", options);
            return;
        }

        // Print modification list
        if(cmd.hasOption("printPTM")){
            ModificationDB.save_mod2file = false;
            ModificationGear.getInstance().printPTM();
            return;
        }

        if(cmd.hasOption("p")){
            String parameter_name = cmd.getOptionValue("p");
            if(parameter_name.equalsIgnoreCase("show")){
                CParameterSet.print_all();
                return;
            }
        }


        if(cmd.hasOption("b")){
            if(cmd.getOptionValue("b").equalsIgnoreCase("show")){
                MSDataSet.print_default_datasets(false);
                return;
            }else if(cmd.getOptionValue("b").equalsIgnoreCase("show_full")){
                MSDataSet.print_default_datasets(true);
                return;
            }
        }

        String out_dir = "." + File.separator;
        if(cmd.hasOption("o")){
            out_dir = cmd.getOptionValue("o");
        }

        if(!cmd.hasOption("db")){
            Cloger.getInstance().logger.error("Protein database (-db) is not provided!");
            System.exit(1);
        }

        if(!cmd.hasOption("i")){
            Cloger.getInstance().logger.error("Input for -i is required!");
            System.exit(1);
        }

        CParameter.cmd = StringUtils.join(args," ");

        if(cmd.hasOption("b")){
            String dataset_id = cmd.getOptionValue("b");
            // use default MS/MS datasets on cloud or public databases
            //CParameter.updateCParameter(cmd);
            HashMap<String,MSDataSet> datasets = MSDataSet.load_datasets(dataset_id);
            if(datasets.size()>=1) {
                Cloger.getInstance().logger.info("The number of MS/MS datasets selected: " + datasets.size() + ", " + String.join(",", datasets.keySet()));
                Collection<MSDataSet> values = datasets.values();
                ArrayList<MSDataSet> datasets_list = new ArrayList<>(values);
                search_multiple_datasets(datasets_list, cmd, out_dir);
            }else{
                CParameter.updateCParameter(cmd);
                // download data from public databases.
                Cloger.getInstance().logger.info("Download data from public database: "+cmd.getOptionValue("b"));
                BuildMSlibrary buildMSlibrary = new BuildMSlibrary();
                if(cmd.hasOption("cpu")){
                    buildMSlibrary.ncpu = Integer.parseInt(cmd.getOptionValue("cpu"));
                }
                String file_type = "";
                if(dataset_id.contains(":")){
                    String[] di = dataset_id.split(":");
                    dataset_id = di[0];
                    file_type = di[1];
                }
                FileDownload.do_raw_data_conversion = true;
                FileDownload.rawDataConvert.delete_raw_file_after_conversion = true;
                FileDownload.rawDataConvert.ms_format = "mgf";
                FileDownload.rawDataConvert.do_gz = false;
                if(cmd.hasOption("dc")){
                    FileDownload.rawDataConvert.convert_bin = cmd.getOptionValue("dc");
                }
                String ms_file = buildMSlibrary.buildLibrary(dataset_id, out_dir, file_type, false);
                search(ms_file,out_dir);
                if(cmd.hasOption("x")){
                    add_extra_validation(cmd.getOptionValue("o")+File.separator+"psm_rank.txt");
                }
            }
        }else if(cmd.hasOption("ms")){
            String ms_file = cmd.getOptionValue("ms");
            CParameter.updateCParameter(cmd);
            search(ms_file,out_dir);
            if(cmd.hasOption("x")){
                add_extra_validation(cmd.getOptionValue("o")+File.separator+"psm_rank.txt");
            }
        }else{
            Cloger.getInstance().logger.warn("Please provide valid parameters!");
        }

    }

    public static void search(String ms_file, String out_dir) throws IOException, SQLException, InterruptedException, ClassNotFoundException {
        CParameter.init();
        SpectraInput.clear();

        String db = CParameter.db;
        File DIR = new File(out_dir);
        if(!DIR.isDirectory()){
            DIR.mkdirs();
        }

        ModificationDB.save_mod2file = true;
        ModificationDB.out_dir = out_dir;
        ModificationDB.getInstance();
        ModificationDB.getInstance().clear_peptide_index();

        CParameter.outdir = DIR.getAbsolutePath()+"/";


        //CParameter.init();
        Cloger.getInstance().logger.info("Start analysis");

        CParameter.print();

        if (CParameter.indexType == 2) {
            Cloger.getInstance().logger.info("Spectrum ID type:" + CParameter.indexType + ", use spectrum title in MS/MS file (MGF) as index for a spectrum.");
        } else {
            Cloger.getInstance().logger.info("Spectrum ID type:" + CParameter.indexType + ", use 1-based number as index for a spectrum.");
        }

        Cloger.getInstance().logger.info("Step 1: target peptide sequence preparation and initial filtering ...");

        String peptideSequence;
        ArrayList<String> pepSeqs = new ArrayList<>();
        if(CParameter.input_target_type.equals(CParameter.TargetInputType.peptide)) {
            ArrayList<CPeptide> cPeptides = new ArrayList<>();
            // take peptide as input
            peptideSequence = CParameter.input_target_file;
            File pFile = new File(peptideSequence);
            if (pFile.isFile()) {
                // Input for "-pep" is a file which contains peptide sequence(s)
                HashMap<String, HashMap<String,TargetedPSM>> pep2targeted_spectra = new HashMap<>();
                HashMap<String, HashMap<String,TargetedPSM>> spectra2targeted_peptide = new HashMap<>();
                HashSet<String> removeDupPep = new HashSet<>();
                BufferedReader pepReader = new BufferedReader(new FileReader(pFile));
                String line;
                boolean score_targeted_psm = false;
                while ((line = pepReader.readLine()) != null) {

                    String[] pep = line.split("\t");
                    String pSeq = pep[0].toUpperCase();
                    if(pep.length>=2){
                        score_targeted_psm = true;
                        // if there are multiple rows in the file, the second column should be spectrum title of interest
                        String spectrum_title = pep[1];
                        if(pep2targeted_spectra.containsKey(pSeq)){
                            pep2targeted_spectra.get(pSeq).put(spectrum_title,new TargetedPSM());
                        }else{
                            HashMap<String, TargetedPSM> spectrumSet = new HashMap<>();
                            spectrumSet.put(spectrum_title, new TargetedPSM());
                            pep2targeted_spectra.put(pSeq,spectrumSet);
                        }

                        if(spectra2targeted_peptide.containsKey(spectrum_title)){
                            spectra2targeted_peptide.get(spectrum_title).put(pSeq,new TargetedPSM());
                        }else{
                            HashMap<String, TargetedPSM> pSet = new HashMap<>();
                            pSet.put(pSeq, new TargetedPSM());
                            spectra2targeted_peptide.put(spectrum_title,pSet);
                        }
                    }
                    if (!removeDupPep.contains(pSeq)) {
                        if(!pep[0].isEmpty()) {
                            // pepSeqs.add(pep[0].toUpperCase());
                            CPeptide cPeptide = new CPeptide();
                            cPeptide.peptideSequence = pSeq;
                            cPeptides.add(cPeptide);
                            removeDupPep.add(pSeq);
                        }
                    }
                }

                pepReader.close();

                if(score_targeted_psm){
                    Cloger.getInstance().logger.info("Spectrum title is provided in file:"+peptideSequence+". Will only score peptide(s) against the spectra provided!");
                    SpectraInput.pep2targeted_spectra = pep2targeted_spectra;
                    SpectraInput.spectra2targeted_peptide = spectra2targeted_peptide;
                    // When provide spectra of interest for comparison, use a low score cutoff.
                    // CParameter.minScore = 4;
                }
            } else {
                String[] peps = peptideSequence.split(",");
                HashSet<String> pepSet = new HashSet<>();
                for(String pep:peps) {
                    if(pepSet.contains(pep)){
                        continue;
                    }else{
                        pepSet.add(pep);
                    }
                    CPeptide cPeptide = new CPeptide();
                    cPeptide.peptideSequence = pep.toUpperCase();
                    // pepSeqs.add(peptideSequence.toUpperCase());
                    cPeptides.add(cPeptide);
                }
            }

            Cloger.getInstance().logger.info("Input peptide sequences:" + cPeptides.size());

            if(CParameter.search_type.equals(CParameter.SearchType.novel)) {
                // if the input peptides are novel peptides
                // if the search type is for known peptide validation, we don't need this step.
                InputProcessor inputProcessor = new InputProcessor();
                inputProcessor.referenceProDB = db;
                inputProcessor.searchRefDB(cPeptides);
                for (CPeptide cPeptide : cPeptides) {
                    if (cPeptide.nhit == 0) {
                        pepSeqs.add(cPeptide.peptideSequence.toUpperCase());
                    } else {
                        Cloger.getInstance().logger.warn("Ignore peptide (reason: exist in reference database): " + cPeptide.peptideSequence);
                    }
                }
            }else{
                for (CPeptide cPeptide : cPeptides) {
                    pepSeqs.add(cPeptide.peptideSequence.toUpperCase());
                }
            }


        }else{
            // take protein, DNA or VCF as input
            String inputSeq = CParameter.input_target_file;
            InputProcessor inputProcessor = new InputProcessor();

            if(CParameter.search_type.equals(CParameter.SearchType.known)){
                // Perform target protein identification. In this case, the value for "-i" is a protein ID.
                if(CParameter.input_target_type.equals(CParameter.TargetInputType.protein)) {
                    Cloger.getInstance().logger.info("Known protein validation ...");
                    TargetProteinID targetProteinID = new TargetProteinID();
                    if (CParameter.add_decoy) {
                        targetProteinID.only_identity_decoy_version = true;
                    }
                    targetProteinID.prepareDB(db, inputSeq, CParameter.outdir);
                    db = targetProteinID.getRefdb();
                    if(targetProteinID.getTargetProteinSequence() ==null || targetProteinID.getTargetProteinSequence().isEmpty()){
                        Cloger.getInstance().logger.error("Target protein is not found in database: "+inputSeq);
                        // use System.exit not return
                        // this is independent to datasets.
                        System.exit(1);
                    }
                    inputSeq = targetProteinID.getTargetProteinSequence();

                }else{
                    Cloger.getInstance().logger.error("We don't support this type of search yet!");
                    // use System.exit not return
                    // this is independent to datasets.
                    System.exit(1);
                }
            }
            inputProcessor.setInputType(CParameter.input_target_type);
            inputProcessor.frames = CParameter.dna_frame;
            inputProcessor.annotationFolder = CParameter.annotation_data_dir;
            inputProcessor.outdir = CParameter.outdir;
            inputProcessor.referenceProDB = db;
            CParameter.internal_db = db;
            ArrayList<CPeptide> cPeptides = inputProcessor.run(inputSeq);
            for(CPeptide cPeptide: cPeptides){
                if(cPeptide.nhit==0) {
                    pepSeqs.add(cPeptide.peptideSequence.toUpperCase());
                }
            }
            int n_shared_peptides = cPeptides.size() - pepSeqs.size();
            Cloger.getInstance().logger.info("Total target peptides:"+cPeptides.size()+", unique peptides:"+pepSeqs.size()+", shared peptides:"+n_shared_peptides);

        }

        if(pepSeqs.size() == 0){
            Cloger.getInstance().logger.info("No valid peptide!");
            // use System.exit not return
            // System.exit(0);
            // use return not System.exit: when searching multiple datasets, different datasets may use different enzyme+miss cleavages preference, this may result in
            // different unique peptides
            return;
        }else{
            Cloger.getInstance().logger.info("Valid target peptides: "+pepSeqs.size());
        }


        ArrayList<PeptideInput> peptideInputs = generatePeptideInputs(pepSeqs);

        Cloger.getInstance().logger.info("Step 1: target peptide sequence preparation and initial filtering done: time elapsed = "+Cloger.getInstance().getRunTime());

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Find candidate spectra for each target peptide
        SpectraInput.out_dir = out_dir;
        Cloger.getInstance().logger.info("Step 2: candidate spectra retrieval and PSM scoring ...");
        File MS_INPUT = new File(ms_file);
        if(MS_INPUT.isDirectory()) {
            Cloger.getInstance().logger.info("Input for MS/MS data is a folder (MS/MS library):" + ms_file);
            SpectraInput.readSpectraFromMSMSlibrary(ms_file, peptideInputs);
        }else{
            //System.out.println(ms_file);
            if(ms_file.startsWith("s3:")){
                if(ms_file.toLowerCase().endsWith(".mgf") ||
                        ms_file.toLowerCase().endsWith(".mgf.gz") ||
                        ms_file.toLowerCase().endsWith(".mgf.xz") ||
                        ms_file.toLowerCase().endsWith(".mzml") ||
                        ms_file.toLowerCase().endsWith(".mzml.gz") ||
                        ms_file.toLowerCase().endsWith(".mzml.xz")
                ){
                    // A file
                    SpectraInput.readMSMSfast(ms_file, peptideInputs);

                }else{
                    // A folder
                    SpectraInput.readSpectraFromMSMSlibrary(ms_file, peptideInputs);
                }
            }else if (ms_file.toLowerCase().endsWith(".mgf")) {
                // A file
                SpectraInput.readMSMSfast(ms_file, peptideInputs);
            }else if(ms_file.startsWith("mzspec:")){
                // USI
                Cloger.getInstance().logger.info("Input MS/MS data is a USI:" + ms_file);
                FileDownload fileDownload = new FileDownload();
                ms_file = fileDownload.download_usi(ms_file,out_dir);
                SpectraInput.readMSMSfast(ms_file, peptideInputs);

            } else {
                // SQL DB
                SpectraInput.readMSMSdb2(ms_file, peptideInputs);
            }
        }

        Cloger.getInstance().logger.info("Time elapsed: "+ Cloger.getInstance().getRunTime());

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate scores for target peptide and the matched spectra if PSMs are not scored.
        if(!SpectraInput.isTargetPeptide2spectrumScored) {
            Cloger.getInstance().logger.info("Start to score PSMs ...");

            ExecutorService fixedThreadPoolScore = Executors.newFixedThreadPool(CParameter.cpu);

            for (PeptideInput peptideInput : peptideInputs) {
                fixedThreadPoolScore.execute(new ScoreWorker(peptideInput));
            }
            fixedThreadPoolScore.shutdown();
            fixedThreadPoolScore.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);
            Cloger.getInstance().logger.info("Start to score PSMs done.");
        }

        Cloger.getInstance().logger.info("Step 2: candidate spectra retrieval and PSM scoring done: time elapsed = "+Cloger.getInstance().getRunTime());

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // search all the spectra against the reference protein database.
        Cloger.getInstance().logger.info("Step 3-4: competitive filtering based on reference sequences and statistical evaluation ...");

        DatabaseInput databaseInput = new DatabaseInput(db);
        databaseInput.add_target_peptides(pepSeqs);
        if(CParameter.use_tag_file){
            databaseInput.tagFile = CParameter.tag_file;
        }

        //if (cmd.hasOption("fixMod")) {
        String[] ifixmod = CParameter.fixMods.split(",");
        for (String i : ifixmod) {
            int id = Integer.parseInt(i);
            // 0 means no modification
            if(id != 0) {
                databaseInput.addFixedModifications(id);
            }
        }
        //}

        //if (cmd.hasOption("varMod")) {
        String[] ivarmod = CParameter.varMods.split(",");
        for (String i : ivarmod) {
            int id = Integer.parseInt(i);
            databaseInput.addVarModifications(id);
        }

        if(db.toLowerCase().endsWith(".fa") || db.toLowerCase().endsWith(".fasta")) {

            String sqldb_file = db+".sqldb";
            File SQLDB = new File(sqldb_file);
            if(SQLDB.isFile()){
                Cloger.getInstance().logger.info("Use indexed database:"+sqldb_file);
                databaseInput.db = sqldb_file;
                databaseInput.readDBSQL();
            }else {
                Cloger.getInstance().logger.info("Don't find indexed database:"+sqldb_file);
                if(CParameter.use_tag_file) {
                    Cloger.getInstance().logger.info("Use tag file:" + CParameter.tag_file);
                    databaseInput.tagFile = CParameter.tag_file;
                    databaseInput.readTag();

                }else{
                    Cloger.getInstance().logger.info("Use database:" + db);
                    //databaseInput.readDB();
                    databaseInput.db_search(peptideInputs);


                }
            }
        }else{
            System.err.println("Please provide valid protein database: .fa or .fasta format.");
            // use System.exit not return.
            // This is independent to datasets.
            System.exit(0);
        }

        Cloger.getInstance().logger.info("Time elapsed: "+Cloger.getInstance().getRunTime());

        Iterator<PeptideInput> inputIterator = peptideInputs.iterator();
        while(inputIterator.hasNext()){
            PeptideInput peptideInput = inputIterator.next();
            boolean rm = true;

            Iterator<JPeptide> jPeptideIterator = peptideInput.getPtmIsoforms().iterator();
            while(jPeptideIterator.hasNext()){
                JPeptide jPeptide = jPeptideIterator.next();
                if(jPeptide.spectraIndexs.size()==0){
                    jPeptideIterator.remove();
                }else{
                    rm = false;
                }
            }
            if(rm){
                inputIterator.remove();
            }
        }

        Cloger.getInstance().logger.info("Peptides with matched spectra: "+peptideInputs.size());

        if(peptideInputs.size()==0){
            Cloger.getInstance().logger.info("Done!");
            // If there is no spectrum found for any peptide, don't use System.exit, instead using return.
            // System.exit(0);
            return;
        }

        // p-value calculation
        StatisticScoreCalc statisticScoreCalc = new StatisticScoreCalc(peptideInputs);
        statisticScoreCalc.run();

        String outfile = export_psm_file(peptideInputs);
        String psm_rank_file = RankPSM.rankPSMs(outfile);


        Cloger.getInstance().logger.info("Step 3-4: competitive filtering based on reference sequences and statistical evaluation done: time elapsed = "+Cloger.getInstance().getRunTime());

        exportMGF(psm_rank_file,peptideInputs);

        if(CParameter.do_modification_search_validation){

            Cloger.getInstance().logger.info("Step 5: competitive filtering based on unrestricted post-translational modification searching ...");
            String ptm_file = ModificationDB.getInstance().doPTMValidation(psm_rank_file, peptideInputs, db, databaseInput.fixedModifications, CParameter.outdir);
            summariseResult(psm_rank_file,ptm_file);
            // annotate best PSM from unrestricted modification searching
            if(CParameter.generate_psm_annotation_data) {
                GeneratePSMAnnotation generatePSMAnnotation = new GeneratePSMAnnotation();
                generatePSMAnnotation.annotateBestPsmFromModSearch(ptm_file);
            }

            Cloger.getInstance().logger.info("Step 5: competitive filtering based on unrestricted post-translational modification searching done: time elapsed = "+Cloger.getInstance().getRunTime());


        }

        if(CParameter.generate_psm_annotation_data) {
            GeneratePSMAnnotation generatePSMAnnotation = new GeneratePSMAnnotation();
            generatePSMAnnotation.run(peptideInputs);
        }


        if(CParameter.plot_spectrum){
            String outfig = CParameter.outdir + "/psm.pdf";
            PeakAnnotation.plot(CParameter.outdir + File.separator + "psm_annotation.txt",outfig,CParameter.intensityThreshold);
        }

        //if(msmsfile.toLowerCase().endsWith(".mgf")){
        //    System.out.println("Add scan number into psm_rank.txt file ...");
        //    addScanNumber2PSM(psm_rank_file);
        //}

        add_delta_score(out_dir);

        if(SpectraInput.pep2targeted_spectra.size()>=1) {
            String type_out_file = CParameter.outdir + File.separator + "psm_type.txt";
            export_targeted_psm_status(psm_rank_file, SpectraInput.pep2targeted_spectra, pepSeqs, type_out_file);
        }

        Cloger.getInstance().logger.info("Time elapsed: "+Cloger.getInstance().getRunTime());
        Cloger.getInstance().logger.info("End.");


    }

    /**
     *
     * @param psm_file
     */
    public static void add_extra_validation(String psm_file) throws IOException, SQLException, InterruptedException, ClassNotFoundException {
        // first get confident identified peptides and generate a new mgf file which only contains confidently
        // identified spectra

        HashSet<String> target_peptides = new HashSet<>();
        HashSet<String> spectrum_titles = new HashSet<>();

        File psmF = new File(psm_file);
        if(!psmF.exists()){
            Cloger.getInstance().logger.info("File doesn't exist:"+psm_file);
            return;
        }

        BufferedReader bReader = new BufferedReader(new FileReader(psmF));

        String head = bReader.readLine();

        head = head.trim();
        String h[] = head.split("\t");
        HashMap<String,Integer> hIndex = new HashMap<>();
        for(int i=0;i<h.length;i++){
            hIndex.put(h[i],i);
        }
        String line;
        while((line = bReader.readLine())!=null){
            line = line.trim();
            String d[] = line.split("\t");

            String conf = d[hIndex.get("confident")];
            if(conf.equalsIgnoreCase("Yes")){
                String msID = d[hIndex.get("spectrum_title")];
                String peptide = d[hIndex.get("peptide")];
                target_peptides.add(peptide);
                spectrum_titles.add(msID);
            }
        }

        bReader.close();

        if(target_peptides.size()>=1){

            String out_dir = psmF.getParent();
            out_dir = out_dir + File.separator + "extra_score";
            File outD = new File(out_dir);
            if(!outD.isDirectory()){
                outD.mkdirs();
            }

            String new_ms_file = out_dir + File.separator + "input.mgf";
            BufferedWriter bWriter = new BufferedWriter(new FileWriter(new_ms_file));

            // get spectra
            String ms_file = psm_file.replaceAll(".txt$",".mgf");
            File mgfFile = new File(ms_file);
            //spectrumFactory.addSpectra(mgfFile, null);
            WaitingHandler waitingHandler = new WaitingHandlerCLIImpl();
            waitingHandler.setDisplayProgress(false);

            MgfFileIterator mgfFileIterator = new MgfFileIterator(mgfFile, waitingHandler);
            String title;
            Spectrum spectrum;
            while ((title = mgfFileIterator.next()) != null) {
                spectrum = mgfFileIterator.getSpectrum();

                if(spectrum_titles.contains(title)){
                    spectrum.spectrumTitle = title;
                    bWriter.write(IndexWorker.asMgf(spectrum,title,spectrum.getPrecursor().possibleCharges[0],null));
                }

            }
            bWriter.close();

            CParameter.TargetInputType save_input_target_type = CParameter.input_target_type;
            String save_input_target_file = CParameter.input_target_file;
            int save_score_method = CParameter.scoreMethod;
            int save_indexType = CParameter.indexType;
            String save_db = CParameter.db;
            CParameter.input_target_type = CParameter.TargetInputType.peptide;
            CParameter.input_target_file = String.join(",",target_peptides);
            CParameter.db = CParameter.internal_db;
            CParameter.indexType = 2;
            if(CParameter.scoreMethod == 1){
                CParameter.scoreMethod = 2;
            }else if(CParameter.scoreMethod == 2){
                CParameter.scoreMethod = 1;
            }


            search(new_ms_file,out_dir);

            CParameter.input_target_type = save_input_target_type;
            CParameter.input_target_file = save_input_target_file;
            CParameter.scoreMethod = save_score_method;
            CParameter.db = save_db;
            CParameter.indexType = save_indexType;
        }


    }
    public static void search_multiple_datasets(ArrayList<MSDataSet> msDataSets, CommandLine cmd_option, String out_dir) {
        if(msDataSets.size()>=1){
            int n_left = msDataSets.size();
            for(int i=0;i<msDataSets.size();i++){
                MSDataSet msDataSet = msDataSets.get(i);
                n_left = n_left - 1;
                Cloger.getInstance().logger.info("Searching MS/MS dataset: " + msDataSet.getDataset_name() + ". "+n_left+" left, "+i+" finished.");
                CParameter.updateCParameter(msDataSet);
                CParameter.updateCParameter(cmd_option);
                String ms_file = msDataSet.getMs_file();
                String out_dir_i = out_dir + "/" + msDataSet.getDataset_name();
                try {
                    search(ms_file,out_dir_i);
                } catch (IOException | SQLException | InterruptedException | ClassNotFoundException e) {
                    e.printStackTrace();
                }

                if(CParameter.add_extra_validation){
                    try {
                        add_extra_validation(out_dir_i+File.separator+"psm_rank.txt");
                    } catch (IOException | SQLException | InterruptedException | ClassNotFoundException e) {
                        e.printStackTrace();
                    }
                }
            }
        }else{
            Cloger.getInstance().logger.warn("No valid MS/MS dataset found!");
        }
    }

    public static void addScanNumber2PSM(String psm_rank_file) throws IOException {

        String old_psm_rank_file = psm_rank_file.replaceAll("txt$","") + "tmp";
        FileUtils.copyFile(new File(psm_rank_file), new File(old_psm_rank_file));

        BufferedReader psm = new BufferedReader(new FileReader(old_psm_rank_file));
        PrintWriter writer = new PrintWriter(psm_rank_file, "UTF-8");
        String hl= psm.readLine().trim();
        String[] hls = hl.split("\t");
        HashMap<String, Integer> hMap = new HashMap<>();
        for (int i = 0; i < hls.length; i++) {
            hMap.put(hls[i], i);
        }

        writer.println(hl+"\tscan_number");
        String line;
        while((line = psm.readLine())!=null){
            line = line.trim();
            String[] d = line.split("\t");
            String title = d[hMap.get("spectrum_title")];
            //String scan_number = SpectraInput.getSpectrum(title).getScanNumber();
            String scan_number = title;
            writer.println(line+"\t"+scan_number);
        }
        psm.close();
        writer.close();

    }


    /**
     *
     * @param psm_rank_file psm_rank.txt
     * @param ptm_file  ptm.txt. This file is from the unrestricted modification searching.
     */
    static void summariseResult(String psm_rank_file, String ptm_file) throws IOException {
        // P-value <= 0.01, rank == 1, n_db == 0
        int identifiedPSMs = 0;

        HashMap<String,ArrayList<ScoreResult>> peptideRes = new HashMap<>();
        HashMap<String,ScoreResult> spectraRes = new HashMap<>();

        String psm_rank_file_header = null;

        try {
            BufferedReader bReader = new BufferedReader(new FileReader(psm_rank_file));

            String line = bReader.readLine().trim();
            psm_rank_file_header = line;
            String[] headLine = line.split("\t");
            HashMap<String,Integer> headMap = new HashMap<>();
            for(int i=0;i<headLine.length;i++){
                headMap.put(headLine[i],i);
            }



            while((line = bReader.readLine())!=null){
                line = line.trim();
                String[] d = line.split("\t");
                String spectrum_title  = d[headMap.get("spectrum_title")];
                String peptideSequence = d[headMap.get("peptide")];

                double pvalue = Double.parseDouble(d[headMap.get("pvalue")]);
                int rank = Integer.parseInt(d[headMap.get("rank")]);
                int n_db = Integer.parseInt(d[headMap.get("n_db")]);
                double score = Double.parseDouble(d[headMap.get("score")]);

                boolean psm_validation_passed = CParameter.passed_psm_validation(peptideSequence.length(),
                        pvalue,n_db,rank);
                if(psm_validation_passed){
                    identifiedPSMs++;
                    ScoreResult scoreResult = new ScoreResult();
                    scoreResult.score = score;
                    scoreResult.peptideSequence = peptideSequence;
                    scoreResult.spectrum_title = spectrum_title;
                    scoreResult.eInfo = line;

                    if(CParameter.scoreMethod == 0){
                        scoreResult.score2 = Double.parseDouble(d[headMap.get("score2")]);
                    }

                    if(peptideRes.containsKey(peptideSequence)){
                        peptideRes.get(peptideSequence).add(scoreResult);
                    }else{
                        ArrayList<ScoreResult> scoreResults = new ArrayList<>();
                        scoreResults.add(scoreResult);
                        peptideRes.put(peptideSequence,scoreResults);
                    }

                    spectraRes.put(spectrum_title,scoreResult);
                }

            }

            bReader.close();

        } catch (IOException e) {
            e.printStackTrace();
        }

        Cloger.getInstance().logger.info("Identified PSMs: "+spectraRes.size());
        Cloger.getInstance().logger.info("Identified peptides: "+peptideRes.size());



        // at least having one PSM with p-value <= 0.01, rank==1 and n_db ==0
        if(identifiedPSMs>=1){
            File PT = new File(ptm_file);
            if(PT.isFile()){
                try {

                    HashMap<String,ArrayList<String>> spectraPtmRes = new HashMap<>();

                    // ptm.txt. This file is from the unrestricted modification searching.
                    BufferedReader bReader = new BufferedReader(new FileReader(PT));

                    String line = bReader.readLine().trim();
                    String[] headLine = line.split("\t");
                    HashMap<String, Integer> headMap = new HashMap<>();
                    for (int i = 0; i < headLine.length; i++) {
                        headMap.put(headLine[i], i);
                    }


                    while ((line = bReader.readLine()) != null) {
                        line = line.trim();
                        String[] d = line.split("\t");
                        String spectrum_title = d[headMap.get("spectrum_title")];
                        //String peptideSequence = d[headMap.get("peptide")];


                        // compare the score from PTM searching with the score from matching to target peptide.

                        if(CParameter.scoreMethod == 0){
                            double score = Double.parseDouble(d[headMap.get("score1")]);
                            double score2 = Double.parseDouble(d[headMap.get("score2")]);

                            // very stringent
                            if(CParameter.unrestrictedSearchWithEqualAndBetterScore) {
                                if (score >= spectraRes.get(spectrum_title).score || score2 >= spectraRes.get(spectrum_title).score2) {

                                    if (spectraPtmRes.containsKey(spectrum_title)) {
                                        spectraPtmRes.get(spectrum_title).add(line);
                                    } else {
                                        ArrayList<String> tmp = new ArrayList<>();
                                        tmp.add(line);
                                        spectraPtmRes.put(spectrum_title, tmp);
                                    }
                                }
                            // it's normal cutoff
                            }else {
                                if (score > spectraRes.get(spectrum_title).score || score2 > spectraRes.get(spectrum_title).score2) {

                                    if (spectraPtmRes.containsKey(spectrum_title)) {
                                        spectraPtmRes.get(spectrum_title).add(line);
                                    } else {
                                        ArrayList<String> tmp = new ArrayList<>();
                                        tmp.add(line);
                                        spectraPtmRes.put(spectrum_title, tmp);
                                    }
                                }
                            }
                        }else{
                            double score = Double.parseDouble(d[headMap.get("score")]);

                            // very stringent
                            if(CParameter.unrestrictedSearchWithEqualAndBetterScore) {
                                if (score >= spectraRes.get(spectrum_title).score) {
                                    if (spectraPtmRes.containsKey(spectrum_title)) {
                                        spectraPtmRes.get(spectrum_title).add(line);
                                    } else {
                                        ArrayList<String> tmp = new ArrayList<>();
                                        tmp.add(line);
                                        spectraPtmRes.put(spectrum_title, tmp);
                                    }
                                }
                            }else {
                                if (score > spectraRes.get(spectrum_title).score) {
                                    if (spectraPtmRes.containsKey(spectrum_title)) {
                                        spectraPtmRes.get(spectrum_title).add(line);
                                    } else {
                                        ArrayList<String> tmp = new ArrayList<>();
                                        tmp.add(line);
                                        spectraPtmRes.put(spectrum_title, tmp);
                                    }
                                }
                            }
                        }

                    }

                    bReader.close();

                    Cloger.getInstance().logger.info("Spectra with higher score from PTM searching: "+spectraPtmRes.size());

                    if(spectraPtmRes.size()>=1) {

                        String ofile = ptm_file;
                        ofile = ofile.replaceAll(".txt$","_detail.txt");
                        BufferedWriter bWriter = new BufferedWriter(new FileWriter(ofile));
                        for (int i = 0; i < headLine.length; i++) {
                            headLine[i] = "ptm_" + headLine[i];
                        }
                        bWriter.write(psm_rank_file_header + "\t" + StringUtils.join(headLine, "\t") + "\n");
                        for (String title : spectraPtmRes.keySet()) {
                            for (String ll : spectraPtmRes.get(title)) {
                                bWriter.write(spectraRes.get(title).eInfo + "\t" + ll+"\n");
                            }
                        }


                        bWriter.close();

                        int nFinalPeps = 0;
                        int nFinalPsms = 0;
                        for(String pep: peptideRes.keySet()){
                            Iterator<ScoreResult> iterator = peptideRes.get(pep).iterator();
                            while(iterator.hasNext()){
                                ScoreResult scoreResult = iterator.next();
                                if(spectraPtmRes.containsKey(scoreResult.spectrum_title)){
                                    iterator.remove();
                                }
                            }

                            if(peptideRes.get(pep).size()>=1){
                                nFinalPeps++;
                                nFinalPsms = nFinalPsms + peptideRes.get(pep).size();
                            }

                        }

                        Cloger.getInstance().logger.info("Identified peptides after PTM searching: "+nFinalPeps);
                        Cloger.getInstance().logger.info("Identified PSMs after PTM searching: "+nFinalPsms);


                    }

                    // Based on the PTM searching result, filtering psm_rank.txt file
                    String old_psm_rank_file = psm_rank_file.replaceAll("txt$","") + "tmp";
                    FileUtils.copyFile(new File(psm_rank_file), new File(old_psm_rank_file));

                    BufferedReader psm = new BufferedReader(new FileReader(old_psm_rank_file));
                    PrintWriter writer = new PrintWriter(psm_rank_file, "UTF-8");
                    String hl= psm.readLine().trim();
                    String[] hls = hl.split("\t");
                    HashMap<String, Integer> hMap = new HashMap<>();
                    for (int i = 0; i < hls.length; i++) {
                        hMap.put(hls[i], i);
                    }
                    // add column confident: Yes/No
                    writer.println(hl+"\tn_ptm\tconfident");

                    while((line = psm.readLine())!=null){
                        line = line.trim();
                        String[] d = line.split("\t");
                        String title = d[hMap.get("spectrum_title")];
                        String peptideSequence = d[hMap.get("peptide")];
                        double pvalue = Double.parseDouble(d[hMap.get("pvalue")]);
                        int rank = Integer.parseInt(d[hMap.get("rank")]);
                        int n_db = Integer.parseInt(d[hMap.get("n_db")]);
                        int n_ptm = -1;
                        boolean psm_validation_passed = CParameter.passed_psm_validation(peptideSequence.length(),
                                pvalue,n_db,rank);

                        // System.out.println(title+"\t"+peptideSequence+"\t"+pvalue+"\t"+psm_validation_passed);

                        if(psm_validation_passed){
                            if(spectraPtmRes.containsKey(title)){
                                n_ptm = spectraPtmRes.get(title).size();
                            }else{
                                n_ptm = 0;
                            }

                        }else{
                            n_ptm = -1;
                        }

                        if(n_ptm==0){
                            writer.println(line+"\t"+n_ptm+"\tYes");
                        }else{
                            writer.println(line+"\t"+n_ptm+"\tNo");
                        }


                    }
                    psm.close();
                    writer.close();


                } catch (IOException e) {
                    e.printStackTrace();
                }

            }else{
                // no unrestricted modification filtering
                String old_psm_rank_file = psm_rank_file.replaceAll("txt$","") + "tmp";
                FileUtils.copyFile(new File(psm_rank_file), new File(old_psm_rank_file));

                BufferedReader psm = new BufferedReader(new FileReader(new File(old_psm_rank_file)));
                PrintWriter writer = new PrintWriter(psm_rank_file, "UTF-8");
                String hl= psm.readLine().trim();
                String hls[] = hl.split("\t");
                HashMap<String, Integer> hMap = new HashMap<>();
                for (int i = 0; i < hls.length; i++) {
                    hMap.put(hls[i], i);
                }
                // add column confident: Yes/No
                writer.println(hl+"\tconfident");
                String line;
                while((line = psm.readLine())!=null){
                    line = line.trim();
                    String d[] = line.split("\t");
                    String title = d[hMap.get("spectrum_title")];
                    double pvalue = Double.parseDouble(d[hMap.get("pvalue")]);
                    int rank = Integer.parseInt(d[hMap.get("rank")]);
                    int n_db = Integer.parseInt(d[hMap.get("n_db")]);

                    String peptideSequence = d[hMap.get("peptide")];
                    boolean psm_validation_passed = CParameter.passed_psm_validation(peptideSequence.length(),
                            pvalue,n_db,rank);

                    if(psm_validation_passed){
                        writer.println(line+"\tYes");
                    }else{
                        writer.println(line+"\tNo");
                    }

                }
                psm.close();
                writer.close();
            }
        }else{
            // no confident identification: the following code is the same with above.
            String old_psm_rank_file = psm_rank_file.replaceAll("txt$","") + "tmp";
            FileUtils.copyFile(new File(psm_rank_file), new File(old_psm_rank_file));

            BufferedReader psm = new BufferedReader(new FileReader(new File(old_psm_rank_file)));
            PrintWriter writer = new PrintWriter(psm_rank_file, "UTF-8");
            String hl= psm.readLine().trim();
            String hls[] = hl.split("\t");
            HashMap<String, Integer> hMap = new HashMap<>();
            for (int i = 0; i < hls.length; i++) {
                hMap.put(hls[i], i);
            }
            if(CParameter.do_modification_search_validation) {
                // add column confident: Yes/No
                if(hMap.containsKey("n_ptm")){
                    writer.println(hl + "\tconfident");
                }else{
                    writer.println(hl + "\tn_ptm\tconfident");
                }

            }else{
                // add column confident: Yes/No
                writer.println(hl + "\tconfident");
            }
            String line;
            while((line = psm.readLine())!=null){
                line = line.trim();
                String[] d = line.split("\t");
                String title = d[hMap.get("spectrum_title")];
                double pvalue = Double.parseDouble(d[hMap.get("pvalue")]);
                int rank = Integer.parseInt(d[hMap.get("rank")]);
                int n_db = Integer.parseInt(d[hMap.get("n_db")]);
                if(CParameter.do_modification_search_validation){
                    writer.println(line + "\t-1\tNo");
                }else {
                    String peptideSequence = d[hMap.get("peptide")];
                    boolean psm_validation_passed = CParameter.passed_psm_validation(peptideSequence.length(),
                            pvalue,n_db,rank);
                    if (psm_validation_passed) {
                        // since there is no confident ID, so this will never happen
                        writer.println(line + "\tYes");
                    } else {
                        writer.println(line + "\tNo");
                    }
                }

            }
            psm.close();
            writer.close();
        }



    }

    public static ScoreResult scorePeptide2Spectrum(Peptide peptide, Spectrum msnSpectrum){
        JSpectrum jSpectrum = new JSpectrum();
        double [] mzArray = msnSpectrum.mz;
        double [] intensityArray = msnSpectrum.intensity;
        for (int i=0;i<mzArray.length;i++) {
            JPeak jPeak = new JPeak(mzArray[i], intensityArray[i]);
            jSpectrum.addRawPeak(jPeak);
        }
        jSpectrum.resetPeaks();
        jSpectrum.sortPeaksByMZ();

        JPeptideSpectrumMatch psm = new JPeptideSpectrumMatch();
        psm.setCalculatedMassToCharge(msnSpectrum.getPrecursor().mz);
        psm.setCharge(msnSpectrum.getPrecursor().possibleCharges[0]);
        psm.setExperimentalMassToCharge(msnSpectrum.getPrecursor().mz);
        psm.setSpectrumID(msnSpectrum.getSpectrumTitle());
        psm.setPepSeq(peptide.getSequence());

        //JMatch.ms2tol = itol;

        double score = -1.0;
        ScoreResult scoreResult = new ScoreResult();

        if(CParameter.scoreMethod == 0){
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // The first scoring algorithm: hyperscore
            HyperscoreMatch hyperscoreMatch = new HyperscoreMatch(msnSpectrum, jSpectrum, psm);
            hyperscoreMatch.objPeptide = peptide;

            try {
                hyperscoreMatch.initialize(false, 2);
            } catch (ClassNotFoundException | IOException | InterruptedException | SQLException e) {
                e.printStackTrace();
            }

            hyperscoreMatch.calcHyperScore();


            score = hyperscoreMatch.getHyperScore();

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // The second scoring algorithm: mvh

            // We need to reset the peaks
            jSpectrum.resetPeaks();
            jSpectrum.sortPeaksByMZ();

            MVHMatch MVHMatch = new MVHMatch(msnSpectrum,jSpectrum, psm);
            MVHMatch.objPeptide = peptide;
            try {
                MVHMatch.initialize(false, 3);
            } catch (ClassNotFoundException | IOException | InterruptedException | SQLException e) {
                e.printStackTrace();
            }
            try {
                MVHMatch.calcMVH();
            } catch (ClassNotFoundException | IOException | InterruptedException | SQLException e) {
                e.printStackTrace();
            }
            double score2 = MVHMatch.getMvh();


            scoreResult.score = score;
            scoreResult.score2 = score2;

        }else if(CParameter.scoreMethod == 1){
            HyperscoreMatch hyperscoreMatch = new HyperscoreMatch(msnSpectrum, jSpectrum, psm);
            hyperscoreMatch.objPeptide = peptide;

            try {
                hyperscoreMatch.initialize(false, 2);
            } catch (ClassNotFoundException | IOException | InterruptedException | SQLException e) {
                e.printStackTrace();
            }

            hyperscoreMatch.calcHyperScore();


            score = hyperscoreMatch.getHyperScore();
            scoreResult.score = score;

        }else if(CParameter.scoreMethod == 2){
            MVHMatch MVHMatch = new MVHMatch(msnSpectrum,jSpectrum, psm);
            MVHMatch.objPeptide = peptide;
            try {
                MVHMatch.initialize(false, 3);
            } catch (ClassNotFoundException | InterruptedException | IOException | SQLException e) {
                e.printStackTrace();
            }
            try {
                MVHMatch.calcMVH();
            } catch (ClassNotFoundException | InterruptedException | IOException | SQLException e) {
                e.printStackTrace();
            }
            score = MVHMatch.getMvh();
            scoreResult.score = score;
        }else{
            System.err.println("Please provide valid score method: 1 or 2. Your input is "+CParameter.scoreMethod);
            // use System.exit not return
            // This is independent to datasets
            System.exit(1);
        }

        //System.out.println("Debug:\t"+peptide.getSequence()+"\t"+psm.getExperimentalMassToCharge()+"\t"+psm.getSpectrumID()+"\t"+String.valueOf(JPeptide.getMass(peptide))+"\t"+ModificationDB.getModificationString(peptide)+"\t"+scoreResult.score);


        return(scoreResult);
    }


    public static void exportMGF(String psm_rank_file, ArrayList<PeptideInput> peptideInputs) throws IOException {

        HashSet<String> spectraTitleSet = new HashSet<>();


        BufferedReader bReader = new BufferedReader(new FileReader(psm_rank_file));

        String line = bReader.readLine().trim();
        String[] headLine = line.split("\t");
        HashMap<String,Integer> headMap = new HashMap<>();
        for(int i=0;i<headLine.length;i++){
            headMap.put(headLine[i],i);
        }

        while((line = bReader.readLine())!=null){
            line = line.trim();
            String[] d = line.split("\t");
            String spectrum_title  = d[headMap.get("spectrum_title")];
            spectraTitleSet.add(spectrum_title);

        }

        bReader.close();

        if(spectraTitleSet.size()>=1){

            HashSet<String> findSpectraSet = new HashSet<>();

            String omgf = psm_rank_file.replaceAll(".txt$",".mgf");

            BufferedWriter mWriter = new BufferedWriter(new FileWriter(omgf));

            for (PeptideInput peptideInput : peptideInputs) {
                //if (peptideInput.getOutLines().size() >= 1) {
                    //String peptideSequence = peptideInput.jPeptide.peptide.getSequence();
                    //for (MSnSpectrum spectrum : peptideInput.jPeptide.mSnSpectrums) {
                    for(JPeptide jPeptide:peptideInput.getPtmIsoforms()){
                        //for (MSnSpectrum spectrum : jPeptide.mSnSpectrums) {
                        for (String title : jPeptide.spectraIndexs) {
                            Spectrum spectrum = SpectraInput.spectraMap.get(title);
                            if (spectraTitleSet.contains(spectrum.getSpectrumTitle())) {
                                if(!findSpectraSet.contains(spectrum.getSpectrumTitle())) {
                                    findSpectraSet.add(spectrum.getSpectrumTitle());
                                    mWriter.write(IndexWorker.asMgf(spectrum,spectrum.getSpectrumTitle(),spectrum.getPrecursor().possibleCharges[0],null));
                                }

                            }
                        }
                    }
                //}
            }

            mWriter.close();

            if(spectraTitleSet.size()==findSpectraSet.size()){
                Cloger.getInstance().logger.info("All spectra have been exported to file: "+omgf);
            }else{
                Cloger.getInstance().logger.warn("Warnings: "+(spectraTitleSet.size()-findSpectraSet.size())+" left!");
            }


        }
    }



    public static void generateSpectraAnnotation(PeptideInput peptideInput, Peptide peptide, Spectrum spectrum){
        String psmMatch_tmp = PeakAnnotation.getPeakAnnotation(peptide, spectrum, true);
        if (psmMatch_tmp != null) {
            peptideInput.addOutPeakAnnotations(psmMatch_tmp);
        }
    }


    public static ArrayList<PeptideInput> generatePeptideInputs(ArrayList<String> pepSeqs) throws InterruptedException {

        Cloger.getInstance().logger.info("Generate peptide objects ...");
        ArrayList<PeptideInput> peptideInputs = new ArrayList<>();

        ModificationGear.getInstance();

        int cpu = CParameter.cpu;

        if(cpu > pepSeqs.size()){
            cpu = pepSeqs.size();
        }

        Cloger.getInstance().logger.info("CPU: "+cpu);
        ExecutorService fixedThreadPoolScore = Executors.newFixedThreadPool(cpu);

        ConcurrentHashMap<String,PeptideInput> peptideInputConcurrentHashMap = new ConcurrentHashMap<>();

        for(String pseq : pepSeqs){
            fixedThreadPoolScore.execute(new CreatePeptideInputWorker(pseq,peptideInputConcurrentHashMap));

        }

        fixedThreadPoolScore.shutdown();

        fixedThreadPoolScore.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);

        for(String pep: peptideInputConcurrentHashMap.keySet()){
            peptideInputs.add(peptideInputConcurrentHashMap.get(pep));
        }

        Cloger.getInstance().logger.info("Generate peptide objects done.");

        return(peptideInputs);
    }


    public static void add_delta_score(String out_dir) throws IOException {
        String psm_rank_file = out_dir + File.separator +"psm_rank.txt";
        File psmF = new File(psm_rank_file);
        if(psmF.exists()){

            String detail_file = out_dir + File.separator +"detail.txt";
            File detailF = new File(detail_file);
            HashMap<String,Double> ref_best_psm = new HashMap<>();
            if(detailF.exists()){
                BufferedReader bReader = new BufferedReader(new FileReader(detailF));
                String head = bReader.readLine();
                head = head.trim();
                String h[] = head.split("\t");
                HashMap<String,Integer> hIndex = new HashMap<>();
                for(int i=0;i<h.length;i++){
                    hIndex.put(h[i],i);
                }
                String line;

                while((line = bReader.readLine())!=null){
                    line = line.trim();
                    String d[] = line.split("\t");
                    String msID = d[hIndex.get("spectrum_title")];
                    double score = Double.parseDouble(d[hIndex.get("score")]);
                    if(ref_best_psm.containsKey(msID)){
                        if(ref_best_psm.get(msID) < score){
                            ref_best_psm.put(msID,score);
                        }
                    }else{
                        ref_best_psm.put(msID,score);
                    }

                }
                bReader.close();
            }

            String mod_detail_file = out_dir + File.separator +"ptm.txt";
            File mod_detailF = new File(mod_detail_file);
            HashMap<String,Double> mod_best_psm = new HashMap<>();
            if(mod_detailF.exists()){
                BufferedReader bReader = new BufferedReader(new FileReader(mod_detailF));
                String head = bReader.readLine();
                head = head.trim();
                String h[] = head.split("\t");
                HashMap<String,Integer> hIndex = new HashMap<>();
                for(int i=0;i<h.length;i++){
                    hIndex.put(h[i],i);
                }
                String line;

                while((line = bReader.readLine())!=null){
                    line = line.trim();
                    String d[] = line.split("\t");
                    String msID = d[hIndex.get("spectrum_title")];
                    double score = Double.parseDouble(d[hIndex.get("score")]);
                    if(mod_best_psm.containsKey(msID)){
                        if(mod_best_psm.get(msID) < score){
                            mod_best_psm.put(msID,score);
                        }
                    }else{
                        mod_best_psm.put(msID,score);
                    }

                }
                bReader.close();
            }


            String old_psm_rank_file = psm_rank_file.replaceAll("txt$","") + "tmp";
            FileUtils.copyFile(new File(psm_rank_file), new File(old_psm_rank_file));

            BufferedReader bReader = new BufferedReader(new FileReader(old_psm_rank_file));
            BufferedWriter bWriter = new BufferedWriter(new FileWriter(psm_rank_file));

            String head = bReader.readLine();

            head = head.trim();
            bWriter.write(head+"\tref_delta_score\tmod_delta_score\n");
            String h[] = head.split("\t");
            HashMap<String,Integer> hIndex = new HashMap<>();
            for(int i=0;i<h.length;i++){
                hIndex.put(h[i],i);
            }
            String line;

            HashMap<String,ArrayList<RankPSM>> psmMap = new HashMap<>();
            while((line = bReader.readLine())!=null){
                line = line.trim();
                String d[] = line.split("\t");

                double score = Double.parseDouble(d[hIndex.get("score")]);
                String msID = d[hIndex.get("spectrum_title")];

                double ref_delta_score = score;
                if(ref_best_psm.containsKey(msID)){
                    ref_delta_score = score - ref_best_psm.get(msID);
                }
                double mod_delta_score = score;
                if(mod_best_psm.containsKey(msID)){
                    mod_delta_score = score - mod_best_psm.get(msID);
                }

                bWriter.write(line+"\t"+ref_delta_score+"\t"+mod_delta_score+"\n");

            }

            bReader.close();
            bWriter.close();
        }

    }


    /**
     * Generate PSM type for each PSM. This is useful when PepQuery is used to validate PSMs.
     * @param psm_rank_file psm_rank.txt file
     * @param out_file output file which contains PSM type information
     * @throws IOException
     */
    public static void export_targeted_psm_status(String psm_rank_file,
                                                  HashMap<String, HashMap<String,TargetedPSM>> pep2targeted_spectra,
                                                  ArrayList<String> valid_peptides,
                                                  String out_file) throws IOException {

        File F = new File(psm_rank_file);
        if(F.isFile()) {

            BufferedReader bReader = new BufferedReader(new FileReader(psm_rank_file));

            String line = bReader.readLine().trim();
            String[] headLine = line.split("\t");
            HashMap<String, Integer> headMap = new HashMap<>();
            for (int i = 0; i < headLine.length; i++) {
                headMap.put(headLine[i], i);
            }

            if(!headMap.containsKey("n_ptm")){
                Cloger.getInstance().logger.warn("n_ptm is not present in "+psm_rank_file);
            }

            // possible situations
            // (1) two different peptides are matched to the same spectrum
            // (2) different peptide forms matched to the same spectrum

            // for each unique peptide sequence and spectrum_title combination,
            // get the highest rank (may not be 1)
            HashMap<String, Integer> pep_spectrum2rank = new HashMap<>();
            while ((line = bReader.readLine()) != null) {
                line = line.trim();
                String[] d = line.split("\t");
                String spectrum_title = d[headMap.get("spectrum_title")];
                String peptide = d[headMap.get("peptide")];
                int rank = Integer.parseInt(d[headMap.get("rank")]);
                String pep_spectrum = peptide + "|" + spectrum_title;
                if (pep_spectrum2rank.containsKey(pep_spectrum)) {
                    if (pep_spectrum2rank.get(pep_spectrum) > rank) {
                        pep_spectrum2rank.put(pep_spectrum, rank);
                    }
                } else {
                    pep_spectrum2rank.put(pep_spectrum, rank);
                }
            }

            bReader.close();


            bReader = new BufferedReader(new FileReader(psm_rank_file));

            // column name
            bReader.readLine();

            HashMap<String,PsmType> psmTypes = new HashMap<>();

            while ((line = bReader.readLine()) != null) {
                line = line.trim();
                String[] d = line.split("\t");
                String spectrum_title = d[headMap.get("spectrum_title")];
                String peptide = d[headMap.get("peptide")];
                int n_db = Integer.parseInt(d[headMap.get("n_db")]);
                double p_value = Double.parseDouble(d[headMap.get("pvalue")]);
                int n_ptm = Integer.parseInt(d[headMap.get("n_ptm")]);
                int rank = Integer.parseInt(d[headMap.get("rank")]);

                PsmType psmType = TargetedPSM.get_psm_status(peptide.length(),rank, n_db, p_value, n_ptm);
                String pep_spectrum = peptide + "|" + spectrum_title;
                if (pep_spectrum2rank.get(pep_spectrum) == rank) {
                    if(psmTypes.containsKey(pep_spectrum)){
                        Cloger.getInstance().logger.error("Duplicated PSM:"+pep_spectrum);
                        // use System.exit not return
                        // if this happens, the analysis should be stopped.
                        System.exit(1);
                    }else{
                        psmTypes.put(pep_spectrum,psmType);
                    }
                }
            }

            bReader.close();
            pep_spectrum2rank.clear();

            BufferedWriter bWriter = new BufferedWriter(new FileWriter(out_file));
            bWriter.write("peptide\tspectrum_title\ttype\n");

            for(String peptide : pep2targeted_spectra.keySet()){
                for(String spectrum_title: pep2targeted_spectra.get(peptide).keySet()){
                    String pep_spectrum = peptide + "|" + spectrum_title;
                    if(psmTypes.containsKey(pep_spectrum)){
                        bWriter.write(peptide+"\t"+spectrum_title+"\t"+psmTypes.get(pep_spectrum).toString()+"\n");
                    }else{
                        if(valid_peptides.contains(peptide)) {
                            if (pep2targeted_spectra.get(peptide).get(spectrum_title).is_status(PsmType.no_spectrum_match)) {
                                bWriter.write(peptide + "\t" + spectrum_title + "\t" + PsmType.no_spectrum_match + "\n");
                            } else {
                                bWriter.write(peptide + "\t" + spectrum_title + "\t" + PsmType.low_score + "\n");
                            }
                        }else{
                            bWriter.write(peptide + "\t" + spectrum_title + "\t" + PsmType.invalid_peptide + "\n");
                        }
                    }
                }
            }
            bWriter.close();
        }

    }

    public static String get_targeted_psm_output(Spectrum spectrum, JPeptide jPeptide, ScoreResult scoreResult){
        StringBuilder outBuilder = new StringBuilder();
        ArrayList<Double> tol_res = CParameter.get_mass_error(spectrum.getPrecursor().getMass(spectrum.getPrecursor().possibleCharges[0]), jPeptide.getMass());
        outBuilder.append(jPeptide.peptide.getSequence())
                .append("\t")
                .append(ModificationDB.getInstance().getModificationString(jPeptide.peptide))
                .append("\t")
                .append(jPeptide.spectraIndexs.size())
                .append("\t")
                .append(spectrum.getSpectrumTitle())
                .append("\t")
                .append(spectrum.getPrecursor().possibleCharges[0])
                .append("\t")
                .append(spectrum.getPrecursor().getMass(spectrum.getPrecursor().possibleCharges[0]))
                .append("\t")
                .append(tol_res.get(0))
                .append("\t")
                .append(tol_res.get(1))
                .append("\t")
                .append(tol_res.get(2))
                .append("\t")
                .append(jPeptide.getMass())
                .append("\t")
                .append(spectrum.getPrecursor().mz)
                .append("\t")
                .append(scoreResult.score)
                .append("\t")
                .append(scoreResult.n_matched_ref_peptides_with_better_score)
                .append("\t")
                .append(scoreResult.n_matched_ref_peptides)
                .append("\t")
                .append(scoreResult.n_random_better_match)
                .append("\t")
                .append(scoreResult.n_total_random_peptides)
                .append("\t")
                .append(scoreResult.p_value);
        return outBuilder.toString();

    }

    public static String export_psm_file(ArrayList<PeptideInput> peptideInputs) throws IOException {
        String outfile = CParameter.outdir + "/psm.txt";
        BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(outfile));
        bufferedWriter.write("peptide\tmodification\tn\tspectrum_title\tcharge\texp_mass\ttol_ppm\ttol_da\tisotope_error\tpep_mass\tmz\tscore\tn_db\ttotal_db\tn_random\ttotal_random\tpvalue\n");
        for(PeptideInput peptideInput: peptideInputs) {
            for (JPeptide jPeptide : peptideInput.getPtmIsoforms()) {
                //for (int i=0;i<jPeptide.mSnSpectrums.size();i++) {
                for (int i = 0; i < jPeptide.spectraIndexs.size(); i++) {
                    String out_line = get_targeted_psm_output(SpectraInput.spectraMap.get(jPeptide.spectraIndexs.get(i)),jPeptide,jPeptide.scores.get(i));
                    bufferedWriter.write(out_line+"\n");
                }
            }
        }
        bufferedWriter.close();
        return outfile;
    }


}
