package main.java.pg;


//import com.amazonaws.auth.AWSCredentials;
import com.amazonaws.auth.AWSStaticCredentialsProvider;
import com.amazonaws.auth.AnonymousAWSCredentials;
import com.amazonaws.auth.BasicAWSCredentials;
//import com.amazonaws.regions.Regions;
import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3ClientBuilder;
import com.amazonaws.services.s3.AmazonS3URI;
import com.amazonaws.services.s3.model.AmazonS3Exception;
import com.amazonaws.services.s3.model.S3Object;
import main.java.util.Cloger;
import org.tukaani.xz.XZInputStream;

import java.io.*;
import java.util.zip.GZIPInputStream;


/**
 * Interact with S3 storage
 */
public class S3Interface {

    private static S3Interface instance = null;

    private AmazonS3 s3Client;

    private S3Interface(){
        //this.s3Client = AmazonS3ClientBuilder.defaultClient();
        s3Client = AmazonS3ClientBuilder
                .standard()
                .withCredentials(new AWSStaticCredentialsProvider(new AnonymousAWSCredentials()))
                .withRegion("us-west-2")
                .build();
        //this.s3Client = AmazonS3ClientBuilder.standard().withRegion(Regions.US_WEST_2).build();
    }

    public static S3Interface getInstance() {
        if (instance == null) {
            instance = new S3Interface();
        }
        return instance;
    }

    public void shutdown(){
        instance.s3Client.shutdown();
    }

    public boolean download(String objPath, String outfile){
        AmazonS3URI uri = new AmazonS3URI(objPath);
        try {
            if (S3Interface.getInstance().getS3Client().doesObjectExist(uri.getBucket(), uri.getKey())) {
                uri = new AmazonS3URI(objPath);
            } else if (S3Interface.getInstance().getS3Client().doesObjectExist(uri.getBucket(), uri.getKey() + ".gz")) {
                objPath = objPath + ".gz";
                uri = new AmazonS3URI(objPath);

            } else if (S3Interface.getInstance().getS3Client().doesObjectExist(uri.getBucket(), uri.getKey() + ".xz")) {
                objPath = objPath + ".xz";
                uri = new AmazonS3URI(objPath);
            }
        }catch (AmazonS3Exception e){
            Cloger.getInstance().logger.error("Error when retrieving MS/MS index file!");
            e.printStackTrace();
            System.exit(1);
        }

        //System.out.println(objPath + " => " + outfile);
        if(S3Interface.getInstance().getS3Client().doesObjectExist(uri.getBucket(),uri.getKey())){
            S3Object s3Object = S3Interface.getInstance().getS3Client().getObject(uri.getBucket(),uri.getKey());
            if(objPath.endsWith(".xz")){
                // xz file
                try {
                    downloadXZip(s3Object.getObjectContent(),outfile);
                } catch (IOException e) {
                    e.printStackTrace();
                }

            }else if(objPath.endsWith(".gz")){
                // gz file
                try {
                    downloadGZip(s3Object.getObjectContent(),outfile);
                } catch (IOException e) {
                    e.printStackTrace();
                }

            }else{
                // txt file
                try {
                    downloadText(s3Object.getObjectContent(),outfile);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            return true;
        }else{
            return false;
        }
    }

    public AmazonS3 getS3Client(){
        return this.s3Client;
    }

    private void downloadGZip(InputStream fileStream, String output) throws IOException {

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

    private void downloadXZip(InputStream fileStream, String output) throws IOException {

        InputStream zipStream = new XZInputStream(fileStream);
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

    private void downloadText(InputStream fileStream, String output) throws IOException {

        Reader decoder = new InputStreamReader(fileStream);
        BufferedReader gReader = new BufferedReader(decoder);

        BufferedWriter pWriter = new BufferedWriter(new FileWriter(new File(output)));
        String line;
        while((line = gReader.readLine())!=null){
            pWriter.write(line+"\n");
        }
        gReader.close();
        pWriter.close();

    }



}
