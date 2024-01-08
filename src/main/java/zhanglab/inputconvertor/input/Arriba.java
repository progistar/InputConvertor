package zhanglab.inputconvertor.input;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.module.ToTranslator;

public class Arriba implements ToTranslator {
	
	public void doTranslate (CommandLine cmd) throws IOException, ParseException {
		Pattern AA_PATTERN = Pattern.compile("[A-Z]*");
		
        String arribaTSVFolder = cmd.getOptionValue("i");
        String inputPattern = cmd.getOptionValue("p");
        String outputPrefix = cmd.getOptionValue("o");
        
        File[] files = new File(arribaTSVFolder).listFiles();
        if(files == null) {
        	System.out.println("-i option should be a path of folder including arriba tsv files");
        	System.exit(1);
        }
        
        BufferedWriter BW = new BufferedWriter(new FileWriter(outputPrefix+".arriba.fasta"));
        int idx = 0;
        for(File file : files) {
        	if(file.getName().startsWith(".")) continue;
        	if(file.getName().endsWith(".tsv") && file.getName().contains(inputPattern)) {
	
	            BufferedReader BR = new BufferedReader(new FileReader(file));
	            String[] header = BR.readLine().split("\t");
	            int peptideIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.ARRIBA_PEPTIDE_FIELD_NAME);
	            String line = null;
	            
	            while((line = BR.readLine()) != null) {
	            	String[] fields = line.split("\t");
	            	String peptide = fields[peptideIdx];
	            	if(peptide.equalsIgnoreCase(".")) continue;
	            	
	            	// FNTLQRRPGWVEYFIAALRGCELVDLADEVASVYQSYQP|rshfvaeagvrwhkrgslqa*
	            	peptide = peptide.toUpperCase();
	            	Matcher matcher = AA_PATTERN.matcher(peptide);
	            	String fullPeptide = "";
	            	while(matcher.find()) {
	            		String partialPeptide = matcher.group();
	            		fullPeptide += partialPeptide;
	            	}
	            	BW.append(">").append(InputConvertorConstants.ARRIBA_HEADER_ID).append("_"+(++idx));
	            	BW.newLine();
	            	BW.append(fullPeptide);
	            	BW.newLine();
	            	
	            }
	            
	            BR.close();
	        }
        }
    	BW.close();
	}
	
}
