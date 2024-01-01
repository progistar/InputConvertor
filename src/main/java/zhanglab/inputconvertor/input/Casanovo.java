package zhanglab.inputconvertor.input;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.data.SimpleMGFSelector;
import zhanglab.inputconvertor.data.Peptide;
import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.module.TopXgInput;

public class Casanovo implements TopXgInput{
	public Casanovo () {}
	///////// Casanovo v3.5.0 index ////////////
	public static int PEPTIDE_INDEX = 1;
	public static int SPECTRA_REF_INDEX = 14;
	public static int SCORE_INDEX = 8;
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
		
		// Options //////////
		
		/**
		 * -i casanovo/
		 * -f mgf/
		 * -p Set01
		 * -o TMT2023_LUAD_Set01
		 * 
		 * For example>
		 * 
		 * target files:
		 * casanovo/TMT_Global_Set01_Fx01.csnv.mztab
		 * casanovo/TMT_Global_Set01_Fx02.csnv.mztab
		 * casanovo/TMT_Global_Set01_Fx03.csnv.mztab
		 * casanovo/TMT_Global_Set01_Fx04.csnv.mztab
		 * mgf/TMT_Global_Set01_Fx01.mgf
		 * mgf/TMT_Global_Set01_Fx02.mgf
		 * mgf/TMT_Global_Set01_Fx03.mgf
		 * mgf/TMT_Global_Set01_Fx04.mgf
		 * 
		 * result files (merging the Fx01~04 files):
		 * casanovo/TMT2023_LUAD_Set01.csnv.ic.mztab
		 */
		/////////////////////
		
        String inputFile = cmd.getOptionValue("i");
        String mgfFileBase = cmd.getOptionValue("f");
        String pattern = cmd.getOptionValue("p");
        String batchId = cmd.getOptionValue("o");
        
        File[] files = new File(inputFile).listFiles();
        if(files == null) {
        	System.out.println("-i option should be a path of folder including .csnv.mztab files");
        	System.exit(1);
        }
        
        BufferedWriter BW = new BufferedWriter(new FileWriter(inputFile+"/"+batchId+".csnv.ic.mztab"));
        String batchHeader = null;
        
        for(File file : files) {
        	if(file.getName().startsWith(".")) continue;
        	if(file.getName().endsWith(".csnv.mztab") && file.getName().contains(pattern)) {
        		System.out.println("read: "+file.getName());
        		SimpleMGFSelector mgf = null;
        		
        		BufferedReader BR = new BufferedReader(new FileReader(file));
        		String line = null;
        		boolean startToRead = false;
        		
        		while((line = BR.readLine()) != null) {
        			// mgf file location
        			if(line.startsWith("MTD") && line.contains("ms_run")) {
        				int len = line.split("\\s")[2].split("\\/").length;
        				String mgfFilePath = line.split("\\s")[2].split("\\/")[len-1];
        				mgfFilePath = mgfFilePath.split("\\.")[0]+".mgf";
        				mgfFilePath = mgfFileBase+"/"+mgfFilePath;
        				
        				File mgfFile = new File(mgfFilePath);
        				mgf = new SimpleMGFSelector(mgfFile);
        			} 
        			// find header
        			else if(line.startsWith("PSH")) {
        				
        		        // building header ///////////////////////////////////////////
        				// if the batch header is already written, then pass
        				if(batchHeader == null) {
        					batchHeader = line;
            				// append header
            				BW.append(batchHeader).append("\t")
            				.append(InputConvertorConstants.IC_SCAN_NUM_FIELD_NAME).append("\t")
            				.append(InputConvertorConstants.IC_TITLE_FIELD_NAME).append("\t")
            				.append(InputConvertorConstants.IC_RT_FIELD_NAME).append("\t")
            				.append(InputConvertorConstants.IC_CHARGE_FIELD_NAME).append("\t")
            				.append(InputConvertorConstants.IC_PEPTIDE_FIELD_NAME);
            				BW.newLine();
        				}
        		        ////////////////////////////////// End of building header ////////////////
        				
        				startToRead = true;
        				
        			} 
        			// convert record
        			else if(startToRead) {
        				String[] fields = line.split("\t");
        				
        				// discard if the score is below than 0
        				if(Double.parseDouble(fields[SCORE_INDEX]) < 0) {
        					continue;
        				}
        				
						// Building record ////////////////////////////////////////
        				
        				Peptide peptide = new Peptide(fields[PEPTIDE_INDEX], InputConvertorConstants.CASANOVO);
        				String spectraRef = fields[SPECTRA_REF_INDEX];
        				int scanIdx = Integer.parseInt(spectraRef.split("\\=")[1]);
        				
        				String title = mgf.scanToTitle.get(scanIdx);
        				String rt = mgf.titleToRT.get(title);
        				
        				int len = title.split("\\.").length;
        				String charge = title.split("\\.")[len-1];
        				
        				
        				for(int i=0; i<fields.length; i++) {
        					BW.append(fields[i]).append("\t");
        				}
        				BW.append(scanIdx+"\t")
        				.append(title).append("\t")
        				.append(rt).append("\t")
        				.append(charge).append("\t")
        				.append(peptide.modPeptide);
        				
        				BW.newLine();
        				
				        /////////////////////////////////// End of building record ///////////////////
        			}
        				
        		}
        		
        		BR.close();
        	}
        }
        BW.close();
		
	}
}
