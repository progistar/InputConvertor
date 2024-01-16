package zhanglab.inputconvertor.input;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import zhanglab.inputconvertor.data.Exon;
import zhanglab.inputconvertor.data.FastaEntry;
import zhanglab.inputconvertor.data.GTFLoader;
import zhanglab.inputconvertor.data.GenomeLoader;
import zhanglab.inputconvertor.data.Transcript;
import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.function.Translator;

public class StringTie {

	public static final int MAX_FLANK_SIZE = 14;
	
	public GTFLoader gtf = null;
	public GTFLoader refGTF = null;
	public GenomeLoader refGenome = null;
	
	public StringTie (File file) throws IOException {
		System.out.println("## StringTie ##");
		System.out.println("Load "+file.getName());
		this.gtf = new GTFLoader(file);
	}
	
	public void enrollGenomeSequence (GenomeLoader refGenome) {
		System.out.println("## StringTie ##");
		System.out.println("Enroll reference genome");
		this.refGenome = refGenome;
	}
	
	private String getTranslation (Transcript transcript, int frame) {
		StringBuilder nucleotide = new StringBuilder();
		StringBuilder protein = new StringBuilder();
		
		
		ArrayList<Exon> exons = transcript.exons;
		
		for(int i=0; i<exons.size(); i++) {
			int eStart = exons.get(i).start - 1;
			int eEnd = exons.get(i).end;
			nucleotide.append(refGenome.getSequence(transcript.chr, eStart, eEnd));
		}

		if(transcript.start.equalsIgnoreCase("+")) {
			protein.append(Translator.translation(nucleotide.toString(), frame));
		} else {
			protein.append(Translator.reverseComplementTranslation(nucleotide.toString(), frame));
		}
		
		return protein.toString();
	}
	
	public ArrayList<FastaEntry> getFastaEntry (double FPKMthreshold) throws IOException {
		ArrayList<FastaEntry> fastaEntries = new ArrayList<FastaEntry>();
		
        this.gtf.geneToTranscripts.forEach((g, ts)->{
    		
    		for(Transcript t : ts) {
    			
    			// fpkm threshold
    			if(t.FPKM <= FPKMthreshold) {
    				continue;
    			}
    			
    			if(t.FPKM == Double.MAX_VALUE) {
    				System.out.println(t.tID);
    			}
    			
    			for(int frame = 0; frame < 3; frame++) {
    				String protein = this.getTranslation(t, frame);
    				// split by X (STOP codon)
    				String[] proteins = protein.split("X");
    				int idx = 1;
    				for(int i=0; i<proteins.length; i++) {
    					protein = proteins[i];
    					
    					// cut
    					if(protein.length() < InputConvertorConstants.MIN_PEPT_LEN) continue;
    					
						FastaEntry entry = new FastaEntry();
						entry.tool = InputConvertorConstants.STRINGTIE_HEADER_ID;
						entry.header = t.tID+"|"+frame+"|"+idx++;
						entry.sequence = protein;
						fastaEntries.add(entry);
    				}
    				
    			}
        	}
        	
        });
        
        return fastaEntries;
	}
}
