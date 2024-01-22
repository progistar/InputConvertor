package zhanglab.inputconvertor.input;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import zhanglab.inputconvertor.data.FastaEntry;
import zhanglab.inputconvertor.data.Transcript;
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
        int gene1Idx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.ARRIBA_GENE1_FIELD_NAME);
        int gene2Idx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.ARRIBA_GENE2_FIELD_NAME);
        int breakPoint1Idx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.ARRIBA_BREAK_POINT1_FIELD_NAME);
        int breakPoint2Idx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.ARRIBA_BREAK_POINT2_FIELD_NAME);
        int frameIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.ARRIBA_FRAME_FIELD_NAME);
        int strand1Idx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.ARRIBA_STRAND1_FIELD_NAME);
        int strand2Idx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.ARRIBA_STRAND2_FIELD_NAME);
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
        	
        	if(fullPeptide.length() < InputConvertorConstants.MIN_PEPT_LEN) continue;
        	
        	String gene1 = fields[gene1Idx];
        	String gene2 = fields[gene2Idx];
        	String bp1 = fields[breakPoint1Idx];
        	String bp2 = fields[breakPoint2Idx];
        	String frame = fields[frameIdx];
        	String strand1 = fields[strand1Idx];
        	String strand2 = fields[strand2Idx];
        	
        	FastaEntry entry = new FastaEntry();
			entry.tool = InputConvertorConstants.ARRIBA_HEADER_ID;
			entry.idx = fastaEntries.size()+1;
			
			//this.tool+this.idx+"|"+this.geneId+"|f:"+this.frame+"|s:"+transcript.strand+"|"+this.description;
			entry.originHeader = gene1+"+"+gene2+"|f:"+frame+"|s:"+strand1+","+strand2+"|"+bp1+","+bp2;
			
			entry.sequence = fullPeptide;
			
			
			entry.frame = 0;
			
			fastaEntries.add(entry);
        	
        }
        
        BR.close();
        
        return fastaEntries;
	}
}
