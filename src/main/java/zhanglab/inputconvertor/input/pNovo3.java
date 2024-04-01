package zhanglab.inputconvertor.input;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.data.SimpleMGFSelector;
import zhanglab.inputconvertor.data.Peptide;
import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.module.TopXgInput;
import zhanglab.inputconvertor.module.TopXgInputGeneric;

public class pNovo3 extends TopXgInputGeneric {
	public pNovo3 () {}
	
	///////// pNovo3 v3.1.5 index ////////////
	public static int CHARGE_IDX = 13;
	public static int PEPTIDE_IDX = 1;
	public static int SCORE_IDX = 2;
	public static String[] FIELDS = {
			"Rank",
			"Peptide",
			"Score",
			"Modification abundace",
			"Precursor mass deviation",
			"Path rank",
			"Score of main ions",
			"Score of internal ions",
			"Continuity score of b ions",
			"Continuity score of y ions",
			"Enzyme score",
			"Std. fragment mass deviations",
			"Max. fragment mass deviations",
			"Charge feature",
			"Raw score",
			"Spearman correlation",
			"Gap feature",
			"Gap feature for N-term"
	};
	////////////////////////////////////////////
	
	/**
	 * -i file path of .csnv.mztab 
	 *    or 
	 *    folder path containing .csnv.mztab files
	 * 
	 * 
	 * @param file.txt
	 * @return
	 * @throws IOException
	 * @throws ParseException 
	 */
	public void topXgInputFormat (String[] args) throws IOException, ParseException {
		parseOptions(args);
        
        File iFile = new File(inputFilePath);
        File sFile = new File(spectrumFilePath);
        File oFile = new File(outputFilePath);
        
        boolean isExsited = oFile.exists();

        BufferedWriter BW = new BufferedWriter(new FileWriter(oFile, isExsited));
        
        // building header ///////////////////////////////////////////
        if(!isExsited) {
	        String batchHeader = "";
	        
	        batchHeader = batchHeader
	        		+InputConvertorConstants.IC_TITLE_FIELD_NAME+"\t"
	        		+InputConvertorConstants.IC_SCAN_NUM_FIELD_NAME+"\t"
	        		+InputConvertorConstants.IC_RT_FIELD_NAME+"\t"
	        		+InputConvertorConstants.IC_CHARGE_FIELD_NAME+"\t"
	        		+InputConvertorConstants.IC_SEARCH_SCORE_FIELD_NAME+"\t"
	        		+InputConvertorConstants.IC_PEPTIDE_FIELD_NAME;
	        
	        for(String header : FIELDS) {
	        	batchHeader += "\t"+header;
	        }
	        
        	BW.append(batchHeader);
            BW.newLine();
        }
        ////////////////////////////////// End of building header ////////////////
        
        System.out.println("read: "+iFile.getName());
		BufferedReader BR = new BufferedReader(new FileReader(iFile));
		SimpleMGFSelector mgf = new SimpleMGFSelector(sFile);
		String line = null;
		
		while((line = BR.readLine()) != null) {
			if(line.startsWith("S")) {
				String[] fields = line.split("\t");
				String title = fields[1].split("\\s")[0];
				
				while((line = BR.readLine()) != null) {
					// each record starts with P[N] (for example, P1, P2, ... P10).
					if(line.startsWith("P")) {
						if(mgf.titleToRT.get(title) == null) {
							continue;
						}
						// Building record ////////////////////////////////////////
						fields = line.split("\t");
						// a is M+oxidation, M+15.995
						Peptide peptide = new Peptide(fields[PEPTIDE_IDX], InputConvertorConstants.PNOVO3);
						// charge feature is different from original MGF. It is supposed to be an error.
						
						int len = title.split("\\.").length;
						String searchScore = fields[SCORE_IDX];
						String scanNum = title.split("\\.")[len-2];
						String charge = title.split("\\.")[len-1];
						
						BW.
						append(title).append("\t").
						append(scanNum).append("\t").
						append(mgf.titleToRT.get(title)).append("\t").
						append(charge).append("\t").
						append(searchScore).append("\t").
						append(peptide.modPeptide);
						
						for(int i=0; i<fields.length; i++) {
							BW.append("\t").append(fields[i]);
						}
						
						BW.newLine();
				        /////////////////////////////////// End of building record ///////////////////
						
					} else {
						break;
					}
				}
			}
		}
		
		BR.close();
        
        BW.close();
	}
	
}
