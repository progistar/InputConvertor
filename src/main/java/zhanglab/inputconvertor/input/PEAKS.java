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

public class PEAKS implements TopXgInput{

	public PEAKS () {}
	
	///////// PEAKS 11 index ////////////
	public static int FILE_IDX = -1;
	public static int SCAN_IDX = -1;
	public static int PEPTIDE_IDX = -1;
	public static int CHARGE_IDX = -1;
	////////////////////////////////////////////
	
	
	public void topXgInputFormat (CommandLine cmd) throws IOException, ParseException {
		/**
		 * -i PEAKS/
		 * -m mgf/
		 * -p Set01
		 * -o TMT2023_LUAD_Set01
		 * 
		 * For example>
		 * 
		 * target files:
		 * PEAKS/TMT_Global_Set01.csv 
		 * mgf/TMT_Global_Set01_Fx01.mgf
		 * mgf/TMT_Global_Set01_Fx02.mgf
		 * mgf/TMT_Global_Set01_Fx03.mgf
		 * mgf/TMT_Global_Set01_Fx04.mgf
		 * 
		 * Only PSMs with "Set01" will be retrived.
		 * 
		 * result files (merging the Fx01~04 files):
		 * pNovo3/TMT2023_LUAD_Set01.PEAKS.ic.tsv
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
        
        BufferedWriter BW = new BufferedWriter(new FileWriter(inputFile+"/"+batchId+".PEAKS.ic.tsv"));
        
        
        String batchHeader = null;
        Hashtable<String, SimpleMGFSelector> mgfFiles = new Hashtable<String, SimpleMGFSelector>();
        
        
        for(File file : files) {
        	if(file.getName().startsWith(".")) continue;
        	if(file.getName().endsWith(".csv") && file.getName().contains(batchPattern)) {
        		System.out.println("read: "+file.getName());
        		BufferedReader BR = new BufferedReader(new FileReader(file));
        		String line = BR.readLine(); // read header
        		
        		// convert to TSV
        		line = replaceCSVtoTSV(line);
        		// building header ///////////////////////////////////////////
                if(batchHeader == null) {
                	batchHeader = line;
                	BW.append(batchHeader).append("\t")
    				.append(InputConvertorConstants.IC_SCAN_NUM_FIELD_NAME).append("\t")
    				.append(InputConvertorConstants.IC_TITLE_FIELD_NAME).append("\t")
    				.append(InputConvertorConstants.IC_RT_FIELD_NAME).append("\t")
    				.append(InputConvertorConstants.IC_CHARGE_FIELD_NAME).append("\t")
    				.append(InputConvertorConstants.IC_PEPTIDE_FIELD_NAME);
    				BW.newLine();
    				
    				String[] header = batchHeader.split("\t");
    				FILE_IDX = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.PEAKS_SOURCE_FILE_FIELD_NAME);
    				PEPTIDE_IDX = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.PEAKS_PEPTIDE_FIELD_NAME);
    				SCAN_IDX = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.PEAKS_SCAN_FIELD_NAME);
    				CHARGE_IDX = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.PEAKS_CHARGE_FIELD_NAME);
                }
                ////////////////////////////////// End of building header ////////////////
        		
        		while((line = BR.readLine()) != null) {
        			line = replaceCSVtoTSV(line);
        			String[] fields = line.split("\t");
        			
        			
        			// Building record ////////////////////////////////////////
        			int lastIdx = fields[FILE_IDX].lastIndexOf(".");
        			String fileName = fields[FILE_IDX].substring(0, lastIdx);
        			
    				String mgfFile = mgfFileBase+"/"+fileName+".mgf";
    				
    				SimpleMGFSelector mgf = mgfFiles.get(mgfFile);
    				if(mgf == null) {
    					mgf = new SimpleMGFSelector(new File(mgfFile));
    					mgfFiles.put(mgfFile, mgf);
    				}
    				// note that if PEAKS runs from .raw files, then the charge state can be altered by PEAKS.
    				String scanNum = fields[SCAN_IDX];
    				String title = mgf.scanToTitle.get(Integer.parseInt(scanNum));
    				String charge = fields[CHARGE_IDX];
        			
    				String rt = mgf.titleToRT.get(title);
    				Peptide peptide = new Peptide(fields[PEPTIDE_IDX], InputConvertorConstants.PEAKS);
    				
    				BW.append(line).append("\t")
    				.append(scanNum).append("\t")
    				.append(title).append("\t")
    				.append(rt).append("\t")
    				.append(charge).append("\t")
    				.append(peptide.modPeptide);
    				
    				BW.newLine();
    				/////////////////////////////////// End of building record ///////////////////
        		}
        		
        		BR.close();
        	}
        }
        
        BW.close();
	}
	
	private String replaceCSVtoTSV (String csv) {
		StringBuilder tsv = new StringBuilder();
		String[] fields = csv.split("\\,");
		
		tsv.append(fields[0].replace("\"", ""));
		for(int i=1; i<fields.length; i++) {
			tsv.append("\t").append(fields[i].replace("\"", ""));
		}
		
		return tsv.toString();
	}
}
