package main.java.db;

import main.java.download.FileDownload;
import main.java.util.Cloger;
import net.sf.jfasta.FASTAElement;
import net.sf.jfasta.FASTAFileReader;
import net.sf.jfasta.impl.FASTAElementIterator;
import net.sf.jfasta.impl.FASTAFileReaderImpl;
import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

public class Database {


    public Database(){

    }


    public static String download_db(String db_source, String species, String version, String out_dir){
        File ODIR = new File(out_dir);
        if(!ODIR.isDirectory()){
            ODIR.mkdirs();
        }

        String out_file = "";
        String db_url = "";
        boolean use_latest_version = version.isEmpty() || version.equalsIgnoreCase("latest") || version.equalsIgnoreCase("current");
        if(db_source.equalsIgnoreCase("Swiss-Prot") || db_source.equalsIgnoreCase("SwissProt")){
            if(species.equalsIgnoreCase("human")){
                if(use_latest_version) {
                    db_url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz";
                    String out_file_gz = out_dir + File.separator + "UP000005640_9606.fasta.gz";
                    out_file = out_dir + File.separator + "UP000005640_9606.fasta";
                    FileDownload.download_single_file(db_url,out_file_gz);
                    try {
                        decompress_gzip(out_file_gz, out_file);
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                }else{
                    Cloger.getInstance().logger.error("Database version is not supported yet:"+version);
                    System.exit(1);
                }
            }

        }else if(db_source.equalsIgnoreCase("UniProt") || db_source.equalsIgnoreCase("UniProtKB")){
            // TODO
        }else if(db_source.equalsIgnoreCase("RefSeq")){
            if(species.equalsIgnoreCase("human")) {
                if (use_latest_version) {
                    db_url = "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_protein.faa.gz";
                    String out_file_gz = out_dir + File.separator + "GRCh38_latest_protein.fasta.gz";
                    FileDownload.download_single_file(db_url,out_file_gz);
                    try {
                        out_file = out_dir + File.separator + "GRCh38_latest_protein.fasta";
                        decompress_gzip(out_file_gz, out_file);
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                } else {
                    Cloger.getInstance().logger.error("Database version is not supported yet:" + version);
                    System.exit(1);
                }
            }
        }else if(db_source.equalsIgnoreCase("gencode")){
            if(species.equalsIgnoreCase("human")) {
                if (use_latest_version) {
                    FileDownload fileDownload = new FileDownload();
                    String ftp_base = "ftp.ebi.ac.uk";
                    String host_extra = "pub/databases/gencode/Gencode_human";
                    String file_pattern = "gencode.*.pc_translations.fa.gz";
                    ArrayList<String> urls = new ArrayList<>();
                    FileDownload.get_file_list_from_ftp(ftp_base,host_extra+"/latest_release",file_pattern,urls);
                    for(String file:urls){
                        if(!file.toLowerCase().contains("grch37_mapping")){
                            // db_url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v40.pc_translations.fa.gz";
                            db_url = file;
                        }
                    }
                    if(db_url.isEmpty()){
                        Cloger.getInstance().logger.error("Doesn't find valid url from ftp://"+ftp_base+"/"+host_extra);
                        System.exit(1);
                    }
                    File F = new File(db_url);
                    String out_file_gz = out_dir + File.separator + F.getName();
                    FileDownload.download_single_file(db_url,out_file_gz);
                    try {
                        out_file = out_dir + File.separator + F.getName().replaceAll("fa.gz$","fasta");
                        decompress_gzip(out_file_gz, out_file);
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                    String format_db = out_file.replaceAll(".fasta$","_format.fasta");
                    try {
                        format_gencode_protein_database(out_file,format_db);
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                    out_file = format_db;

                } else {
                    Cloger.getInstance().logger.error("Database version is not supported yet:" + version);
                    System.exit(1);
                }
            }
        }
        return out_file;

    }

    private static void decompress_gzip(String gzip_file, String output) throws IOException {
        FileInputStream fileStream = new FileInputStream(gzip_file);
        InputStream zipStream = new GZIPInputStream(fileStream);
        Reader decoder = new InputStreamReader(zipStream);
        BufferedReader gReader = new BufferedReader(decoder);

        BufferedWriter pWriter = new BufferedWriter(new FileWriter(new File(output)));
        String line;
        while((line = gReader.readLine())!=null){
            pWriter.write(line+"\n");
        }
        gReader.close();
        pWriter.close();
    }


    private static void format_gencode_protein_database(String original_db, String out_db) throws IOException {
        File dbFile = new File(original_db);
        FASTAFileReader reader = new FASTAFileReaderImpl(dbFile);
        FASTAElementIterator it = reader.getIterator();

        BufferedWriter bWriter = new BufferedWriter(new FileWriter(out_db));
        while (it.hasNext()) {
            FASTAElement el = it.next();
            el.setLineLength(1);
            String headLine[] = el.getHeader().split("\\s+");
            // >ENSP00000317992.6|ENST00000327044.7|ENSG00000188976.11|OTTHUMG000000407
            String proID = headLine[0].split("\\|")[0];
            // System.out.println(proID);
            String proSeq = el.getSequence().toUpperCase();
            bWriter.write(">"+proID + " " + el.getHeader() +"\n"+proSeq+"\n");

        }
        reader.close();
        bWriter.close();
    }



}
