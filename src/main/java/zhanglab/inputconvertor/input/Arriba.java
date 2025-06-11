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
        int gene1Idx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.ARRIBA_GENE1_FIELD_NAME);
        int gene2Idx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.ARRIBA_GENE2_FIELD_NAME);
        int geneId1Idx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.ARRIBA_GENEID1_FIELD_NAME);
        int geneId2Idx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.ARRIBA_GENEID2_FIELD_NAME);
        int breakPoint1Idx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.ARRIBA_BREAK_POINT1_FIELD_NAME);
        int breakPoint2Idx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.ARRIBA_BREAK_POINT2_FIELD_NAME);
        int frameIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.ARRIBA_FRAME_FIELD_NAME);
        int strand1Idx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.ARRIBA_STRAND1_FIELD_NAME);
        int strand2Idx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.ARRIBA_STRAND2_FIELD_NAME);
        int transcript1Idx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.ARRIBA_TRANSCRIPT1_FIELD_NAME);
        int transcript2Idx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.ARRIBA_TRANSCRIPT2_FIELD_NAME);
        
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
        	String geneId1 = fields[geneId1Idx];
        	String geneId2 = fields[geneId2Idx];
        	String bp1 = fields[breakPoint1Idx];
        	String bp2 = fields[breakPoint2Idx];
        	String frame = fields[frameIdx];
        	String strand1 = fields[strand1Idx];
        	String strand2 = fields[strand2Idx];
        	String transcript1 = fields[transcript1Idx];
        	String transcript2 = fields[transcript2Idx];
        	
        	FastaEntry entry = new FastaEntry();
			entry.tool = InputConvertorConstants.ARRIBA_HEADER_ID;
			entry.idx = fastaEntries.size()+1;
			entry.geneName = gene1+"+"+gene2;
			entry.geneId = geneId1+"+"+geneId2;
			entry.transcriptId = transcript1+"+"+transcript2;
			entry.frame = frame;
			entry.strand = strand1+","+strand2;
			entry.description = bp1+","+bp2;
			//entry.originHeader = gene1+"+"+gene2+"|"+transcript1+"+"+transcript2+"|"+gene1+"+"+gene2+"|"+frame+"|"+strand1+","+strand2+"|"+bp1+","+bp2;
			
			entry.sequence = fullPeptide;
			fastaEntries.add(entry);
        	
        }
        
        BR.close();
        
        return fastaEntries;
	}
}
