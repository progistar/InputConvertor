package zhanglab.inputconvertor.input;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import zhanglab.inputconvertor.data.FastaEntry;
import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.function.Translator;

public class MANTA {

	public File MANTAFile = null;
	
	public MANTA (File arribaFile) {
		System.out.println("## Arriba ##");
		System.out.println("Load "+arribaFile.getName());
		this.MANTAFile = arribaFile;
	}
	
	public ArrayList<FastaEntry> getFastaEntry () throws IOException {
		ArrayList<FastaEntry> fastaEntries = new ArrayList<FastaEntry>();
		Pattern contigPattern = Pattern.compile("(CONTIG=[A|C|T|G|N]*)");
		Pattern svPattern = Pattern.compile("(SVTYPE=[A-Z]*)");
		BufferedReader BR = null;
		
		if(this.MANTAFile.getName().endsWith(".gz")) {
			FileInputStream fis = new FileInputStream(this.MANTAFile);
            GZIPInputStream gis = new GZIPInputStream(fis);
            InputStreamReader isr = new InputStreamReader(gis);
			BR = new BufferedReader(isr);
		} else {
			BR = new BufferedReader(new FileReader(this.MANTAFile));
		}
		
		// skip meta
		
        String line = null;
        Hashtable<String, String> contigHashtable = new Hashtable<String, String>();
        
        while((line = BR.readLine()) != null) {
        	// skip meta
        	if(line.startsWith("#")) continue;
        	
        	Matcher matcher = contigPattern.matcher(line.toUpperCase());
        	String contig = null;
        	while(matcher.find()) {
        		contig = matcher.group();
        		contig = contig.split("\\=")[1];
        	}
        	
        	// skip if the contig contains "N"
        	if(contig == null || contig.contains("N")) continue;
        	
        	matcher = svPattern.matcher(line.toUpperCase());
        	String svType = null;
        	while(matcher.find()) {
        		svType = matcher.group();
        		svType = svType.split("\\=")[1];
        	}
        	
        	// skip INDELs
        	if(svType.equalsIgnoreCase("INS") || svType.equalsIgnoreCase("DEL")) {
        		continue;
        	}
        	
        	if(contigHashtable.get(contig) != null) {
        		continue;
        	}
        	
        	contigHashtable.put(contig, "");
        	
        	String[] fields = line.split("\t");
        	String mantaId = fields[2];
        	
        	// do six-frame translation
        	for(int i=0; i<3; i++) {
        		FastaEntry entry = new FastaEntry();
    			entry.tool = InputConvertorConstants.MANTA_HEADER_ID;
    			entry.idx = fastaEntries.size()+1;
    			
    			entry.sequence = Translator.translation(contig, i)[1];
    			entry.strand = "+";
    			entry.frame = i+"";
    			entry.transcriptId = mantaId;
    			entry.geneId = mantaId;
    			entry.geneName = mantaId.split("\\:")[0];
    			entry.description = "SV";
    			fastaEntries.add(entry);
    			
    			entry = new FastaEntry();
    			entry.tool = InputConvertorConstants.MANTA_HEADER_ID;
    			entry.idx = fastaEntries.size()+1;
    			
    			entry.sequence = Translator.reverseComplementTranslation(contig, i)[1];
    			entry.strand = "-";
    			entry.frame = i+"";
    			entry.transcriptId = mantaId;
    			entry.geneId = mantaId;
    			entry.geneName = mantaId.split("\\:")[0];
    			entry.description = "SV";
    			fastaEntries.add(entry);
        	}
        	
        }
        
        BR.close();
        
        return fastaEntries;
	}
}
