package zhanglab.inputconvertor.input;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.data.SimpleMGFSelector;
import zhanglab.inputconvertor.data.Peptide;
import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.module.TopXgInput;

public class pNovo3 implements TopXgInput{
	public pNovo3 () {}
	
	///////// pNovo3 v3.1.5 index ////////////
	public static int CHARGE_IDX = 13;
	public static int PEPTIDE_IDX = 1;
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
	public void topXgInputFormat (CommandLine cmd) throws IOException, ParseException {
		/**
		 * -i pNovo3/
		 * -m mgf/
		 * -p Set01
		 * -o TMT2023_LUAD_Set01
		 * 
		 * For example>
		 * 
		 * target files:
		 * pNovo3/pNovo.res (** the file contains PSMs from TMT_Global_Set01 or TMT_Global_Set02)
		 * mgf/TMT_Global_Set01_Fx01.mgf
		 * mgf/TMT_Global_Set01_Fx02.mgf
		 * mgf/TMT_Global_Set01_Fx03.mgf
		 * mgf/TMT_Global_Set01_Fx04.mgf
		 * 
		 * Only PSMs with "Set01" will be retrived.
		 * 
		 * result files (merging the Fx01~04 files):
		 * pNovo3/TMT2023_LUAD_Set01.pNovo3.ic.tsv
		 */
		/////////////////////
        String inputFile = cmd.getOptionValue("i");
        String mgfFileBase = cmd.getOptionValue("f");
        String batchPattern = cmd.getOptionValue("p");
        String batchId = cmd.getOptionValue("o");
        
        File[] files = new File(inputFile).listFiles();
        if(files == null) {
        	System.out.println("-i option should be a path of folder including .res files");
        	System.exit(1);
        }
        
        BufferedWriter BW = new BufferedWriter(new FileWriter(inputFile+"/"+batchId+".pNovo3.ic.tsv"));
        
        // building header ///////////////////////////////////////////
        String batchHeader = "";
        for(String header : FIELDS) {
        	batchHeader += header+"\t";
        }
        batchHeader+=InputConvertorConstants.IC_SCAN_NUM_FIELD_NAME+"\t"
        		+InputConvertorConstants.IC_TITLE_FIELD_NAME+"\t"
        		+InputConvertorConstants.IC_RT_FIELD_NAME+"\t"
        		+InputConvertorConstants.IC_CHARGE_FIELD_NAME+"\t"
        		+InputConvertorConstants.IC_PEPTIDE_FIELD_NAME;
        BW.append(batchHeader);
        BW.newLine();
        ////////////////////////////////// End of building header ////////////////
        
        Hashtable<String, SimpleMGFSelector> mgfFiles = new Hashtable<String, SimpleMGFSelector>();
        for(File file : files) {
        	if(file.getName().startsWith(".")) continue;
        	if(file.getName().endsWith(".res")) {
        		System.out.println("read: "+file.getName());
        		BufferedReader BR = new BufferedReader(new FileReader(file));
        		String line = null;
        		
        		while((line = BR.readLine()) != null) {
        			String[] fields = line.split("\t");
        			if(line.startsWith("S") && line.contains(batchPattern)) {
        				String title = fields[1];
        				String mgfFile = mgfFileBase+"/"+title.split("\\.")[0]+".mgf";
        				
        				SimpleMGFSelector mgf = mgfFiles.get(mgfFile);
        				if(mgf == null) {
        					mgf = new SimpleMGFSelector(new File(mgfFile));
        					mgfFiles.put(mgfFile, mgf);
        				}
        				
        				
        				while((line = BR.readLine()) != null) {
        					// each record starts with P[N] (for example, P1, P2, ... P10).
        					if(line.startsWith("P")) {
        						
        						// Building record ////////////////////////////////////////
        						fields = line.split("\t");
        						// a is M+oxidation, M+15.995
        						Peptide peptide = new Peptide(fields[PEPTIDE_IDX], InputConvertorConstants.PNOVO3);
        						// charge feature is different from original MGF. It is supposed to be an error.
        						
        						int len = title.split("\\.").length;
        						String scanNum = title.split("\\.")[len-2];
        						String charge = title.split("\\.")[len-1];
        						
        						for(int i=0; i<fields.length; i++) {
        							BW.append(fields[i]).append("\t");
        						}
        						BW.append(scanNum).append("\t").
        						append(title).append("\t").
        						append(mgf.titleToRT.get(title)).append("\t").
        						append(charge).append("\t").
        						append(peptide.modPeptide);
        						BW.newLine();
        				        /////////////////////////////////// End of building record ///////////////////
        						
        					} else {
        						break;
        					}
        				}
        			}
        		}
        		
        		BR.close();
        	}
        }
        
        BW.close();
		
	}
}
