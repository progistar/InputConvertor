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

public class Reference {

	public GTFLoader refGTF = null;
	public GenomeLoader refGenome = null;
	
	public Reference (GTFLoader refGTF) throws IOException {
		System.out.println("## Reference ##");
		this.refGTF = refGTF;
	}
	
	public void enrollGenomeSequence (GenomeLoader refGenome) {
		System.out.println("## Reference ##");
		System.out.println("Enroll reference genome");
		this.refGenome = refGenome;
	}
	
	private ArrayList<FastaEntry> getTranslation (Transcript transcript) {
		ArrayList<Exon> exons = transcript.cdss;
		
		refGenome.setSequence(transcript.chr, exons);
		
		ArrayList<FastaEntry> entries = FastaEntry.enumerateFastaEntry(transcript, exons, true);
		
		return entries;
	}
	
	public ArrayList<FastaEntry> getFastaEntry () throws IOException {
		ArrayList<FastaEntry> fastaEntries = new ArrayList<FastaEntry>();
		int[] proteinCodingTranscripts = new int[1];
        this.refGTF.geneToTranscripts.forEach((g, ts)->{
    		
    		for(Transcript t : ts) {
    			if(t.isProteinCoding) {
    				proteinCodingTranscripts[0]++;
    				ArrayList<FastaEntry> entries = this.getTranslation(t);
    				// put gene ID
    				for(FastaEntry entry : entries) {
    					entry.geneId = g;
    					entry.transcriptIds.add(entry.transcript.tID);
    					entry.tool = InputConvertorConstants.REF_HEADER_ID;
    				}
    				fastaEntries.addAll(entries);
    			}
        	}
        	
        });
        
        // remove duplications
        System.out.println("A total of reference protein coding transcripts: "+proteinCodingTranscripts[0]);
        return FastaEntry.removeDuplications(fastaEntries);
	}
}
