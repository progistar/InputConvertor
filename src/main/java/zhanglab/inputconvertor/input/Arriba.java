package zhanglab.inputconvertor.input;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import zhanglab.inputconvertor.data.FastaEntry;
import zhanglab.inputconvertor.env.InputConvertorConstants;

public class Arriba {
	
	public File arribaFile = null;
	
	public Arriba (File arribaFile) {
		System.out.println("## Arriba ##");
		System.out.println("Load "+arribaFile.getName());
		this.arribaFile = arribaFile;
	}
	
	public ArrayList<FastaEntry> getFastaEntry () throws IOException {
		ArrayList<FastaEntry> fastaEntries = new ArrayList<FastaEntry>();
		
		Pattern AA_PATTERN = Pattern.compile("[A-Z]*");
		
		BufferedReader BR = new BufferedReader(new FileReader(this.arribaFile));
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
        	
        }
        
        BR.close();
        
        return fastaEntries;
	}
}
