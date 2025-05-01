package zhanglab.inputconvertor.input;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.data.Peptide;
import zhanglab.inputconvertor.data.SimpleSpectraSelector;
import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.module.TopXgInputGeneric;

public class PiPrimeNovo extends TopXgInputGeneric {
	public PiPrimeNovo () {}
	///////// Casanovo v3.5.0 index ////////////
	public static int PEPTIDE_INDEX = 1;
	public static int SPECTRA_REF_INDEX = 0;
	public static int SCORE_INDEX = 2;
	////////////////////////////////////////////
	/**
	 * 
	 * 
	 * @param file.txt
	 * @return
	 * @throws IOException
	 * @throws ParseException 
	 */
	public void topXgInputFormat (String[] args) throws IOException, ParseException {
		parseOptions(args);
		Hashtable<String, String> ptmTable = loadPTMTable();
        
        File iFile = new File(inputFilePath);
        File sFile = new File(spectrumFilePath);
        File oFile = new File(outputFilePath);
        
        boolean isExsited = oFile.exists();
        
        BufferedWriter BW = new BufferedWriter(new FileWriter(oFile, isExsited));
        String batchHeader = null;
        
        System.out.println("read: "+iFile.getName());
		SimpleSpectraSelector mgf = new SimpleSpectraSelector(sFile);
		BufferedReader BR = new BufferedReader(new FileReader(iFile));
		String line = null;
		boolean startToRead = false;
		
		while((line = BR.readLine()) != null) {
			// find header
			if(!startToRead) {
				
		        // building header ///////////////////////////////////////////
				// if the batch header is already written, then pass
				batchHeader = line;
				// append header
				if(!isExsited) {
					BW
					.append(InputConvertorConstants.IC_TITLE_FIELD_NAME).append("\t")
					.append(InputConvertorConstants.IC_SCAN_NUM_FIELD_NAME).append("\t")
					.append(InputConvertorConstants.IC_RT_FIELD_NAME).append("\t")
					.append(InputConvertorConstants.IC_CHARGE_FIELD_NAME).append("\t")
					.append(InputConvertorConstants.IC_SEARCH_SCORE_FIELD_NAME).append("\t")
					.append(InputConvertorConstants.IC_PEPTIDE_FIELD_NAME).append("\t")
					.append(batchHeader);
					BW.newLine();
				}
		        ////////////////////////////////// End of building header ////////////////
				
				startToRead = true;
				
				String[] headerSplit = batchHeader.split("\t");
				SCORE_INDEX = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.piPRIME_NOVO_SCORE_FIELD_NAME);
				PEPTIDE_INDEX = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.piPRIME_NOVO_PEPTIDE_FIELD_NAME);
				SPECTRA_REF_INDEX = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.piPRIME_NOVO_SPECTRA_REF_FIELD_NAME);
			} 
			// convert record
			else if(startToRead) {
				String[] fields = line.split("\t");
				
				// discard if the score is below than 0
				if(Double.parseDouble(fields[SCORE_INDEX]) < 0) {
					continue;
				}
				
				// Building record ////////////////////////////////////////
				// skip null peptide
				if(fields[PEPTIDE_INDEX].length() == 0) {
					continue;
				}
				Peptide peptide = new Peptide(fields[PEPTIDE_INDEX], ptmTable);
				if(!peptide.isPass) {
					System.out.println("Unknown characters were detected: "+peptide.modPeptide);
					continue;
				}
				
				String spectraRef = fields[SPECTRA_REF_INDEX];
				
				String[] split = spectraRef.split("\\.");
				int scanIdx = Integer.parseInt(split[split.length-2]);
				String charge = split[split.length-1];
				String title = mgf.scanToTitle.get(scanIdx);
				String rt = mgf.titleToRT.get(title);
				String searchScore = fields[SCORE_INDEX];
				
				
				BW
				.append(title).append("\t")
				.append(scanIdx+"\t")
				.append(rt).append("\t")
				.append(charge).append("\t")
				.append(searchScore).append("\t")
				.append(peptide.modPeptide);
				
				for(int i=0; i<fields.length; i++) {
					BW.append("\t").append(fields[i]);
				}
				
				BW.newLine();
				
		        /////////////////////////////////// End of building record ///////////////////
			}
				
		}
		
		BR.close();
        BW.close();
		
	}
	
}
