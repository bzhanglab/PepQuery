package main.java.PSMMatch;
import main.java.util.Cloger;
import java.io.*;
import java.util.Scanner;


/**
 * replace the title names with the index in mgf file 
 *
 */
public class FormatMGF {



	public static void main(String[] args) throws IOException {

		String inputmgf = args[0];
		String outputmgf = args[1];
		BufferedWriter bWriter = new BufferedWriter(new FileWriter(new File(outputmgf)));
		Scanner mgfScanner =new Scanner(new File(inputmgf));			
		mgfScanner.useDelimiter("BEGIN IONS"); 
		int index =-1;
		while(mgfScanner.hasNext()){
			String spectra =mgfScanner.next();
			if(!spectra.contains("END IONS")){
				continue;
			}
			index++;
			if(spectra.contains("TITLE=")){
				spectra=spectra.replaceFirst("TITLE=.*", "TITLE="+index);
				spectra="BEGIN IONS"+spectra;
			}else{
				spectra="BEGIN IONS\n"+"TITLE="+index+"\n"+spectra;
			}
			bWriter.write(spectra);
			
		}
		mgfScanner.close();
		bWriter.close();
	}
	
	
	public static void formatMGF(String inputmgf, String outputmgf) throws IOException{
		Cloger.getInstance().logger.info("Format file:"+inputmgf+" ==> "+outputmgf);
		BufferedWriter bWriter = new BufferedWriter(new FileWriter(new File(outputmgf)));
		Scanner mgfScanner =new Scanner(new File(inputmgf));			
		mgfScanner.useDelimiter("BEGIN IONS"); 
		int index =-1;
		while(mgfScanner.hasNext()){
			String spectra =mgfScanner.next();
			if(!spectra.contains("END IONS")){
				continue;
			}
			index++;
			if(spectra.contains("TITLE=")){
				spectra=spectra.replaceFirst("TITLE=.*", "TITLE="+index);
				spectra="BEGIN IONS"+spectra;
			}else{
				spectra="BEGIN IONS\n"+"TITLE="+index+"\n"+spectra;
			}
			bWriter.write(spectra);
			
		}
		mgfScanner.close();
		bWriter.close();
	}

	
	

}
