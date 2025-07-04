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

public class PEAKS12 extends TopXgInputGeneric {

	public PEAKS12 () {}
	
	///////// PEAKS 11 index ////////////
	public static int FILE_IDX = -1;
	public static int SCAN_IDX = -1;
	public static int PEPTIDE_IDX = -1;
	public static int CHARGE_IDX = -1;
	public static int SCORE_IDX = -1;
	public static int PRECURSOR_IDX = -1;
	public static int MZ_IDX = -1;
	////////////////////////////////////////////
	
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
		BufferedReader BR = new BufferedReader(new FileReader(iFile));
		String line = BR.readLine(); // read header
		
		// convert to TSV
		line = replaceCSVtoTSV(line);
		batchHeader = line;
		// building header ///////////////////////////////////////////
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
		
		String[] header = batchHeader.split("\t");
		FILE_IDX = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.PEAKS_SOURCE_FILE_FIELD_NAME);
		PEPTIDE_IDX = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.PEAKS_PEPTIDE_FIELD_NAME);
		SCAN_IDX = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.PEAKS_SCAN_FIELD_NAME);
		CHARGE_IDX = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.PEAKS_CHARGE_FIELD_NAME);
		
		// first check deep novo score
		SCORE_IDX = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.PEAKS_DEEPNOVO_SCORE_FIELD_NAME);
		if(SCORE_IDX == -1) {
			SCORE_IDX = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.PEAKS_SCORE_FIELD_NAME);
		}
		PRECURSOR_IDX = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.PEAKS_PRECURSOR_ID_NAME);
		MZ_IDX = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.PEAKS_PEPMASS_NAME);
        ////////////////////////////////// End of building header ////////////////
		SimpleSpectraSelector mgf = new SimpleSpectraSelector(sFile);
		String mgfFileName = sFile.getName().substring(0, sFile.getName().lastIndexOf("."));
		
		while((line = BR.readLine()) != null) {
			line = replaceCSVtoTSV(line);
			String[] fields = line.split("\t");
			
			String fileName = fields[FILE_IDX].substring(0, fields[FILE_IDX].lastIndexOf("."));
			if(!mgfFileName.equalsIgnoreCase(fileName)) continue;
			
			
			// Building record ////////////////////////////////////////
			// note that if PEAKS runs from .raw files, then the charge state can be altered by PEAKS.
			String scanNum = fields[SCAN_IDX];
			int mass = (int) (1000 * Double.parseDouble(fields[MZ_IDX]));
			String title = fields[PRECURSOR_IDX] + "+" + mass;
			String charge = fields[CHARGE_IDX];
			String searchScore = fields[SCORE_IDX];
			
			String rt = mgf.titleToRT.get(title);
			Peptide peptide = new Peptide(fields[PEPTIDE_IDX], ptmTable);
			if(!peptide.isPass) {
				System.out.println("Unknown characters were detected: "+peptide.modPeptide);
				continue;
			}
			
			BW
			.append(title).append("\t")
			.append(scanNum).append("\t")
			.append(rt).append("\t")
			.append(charge).append("\t")
			.append(searchScore).append("\t")
			.append(peptide.modPeptide).append("\t")
			.append(line);
			
			BW.newLine();
			/////////////////////////////////// End of building record ///////////////////
		}
		
		BR.close();
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
