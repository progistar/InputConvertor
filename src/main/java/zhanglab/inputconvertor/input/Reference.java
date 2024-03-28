package zhanglab.inputconvertor.input;

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
	
	private ArrayList<FastaEntry> getTranslation (Transcript transcript, boolean isOnlyPC) {
		ArrayList<Exon> exons = null;
		ArrayList<FastaEntry> entries = null;
		if(isOnlyPC) {
			exons = transcript.cdss;
			entries = FastaEntry.enumerateFastaEntryCDS(refGenome, transcript, exons);
		} else {
			exons = transcript.exons;
			entries = FastaEntry.enumerateFastaEntry(refGenome, transcript, exons);
		}
		
		return entries;
	}
	
	public ArrayList<FastaEntry> getFastaEntry (boolean isOnlyPC) throws IOException {
		ArrayList<FastaEntry> fastaEntries = new ArrayList<FastaEntry>();
        this.refGTF.geneToTranscripts.forEach((g, ts)->{
    		
    		for(Transcript t : ts) {
    			if(isOnlyPC && t.isProteinCoding) {
    				ArrayList<FastaEntry> entries = this.getTranslation(t, isOnlyPC);
    				// put gene ID
    				for(FastaEntry entry : entries) {
    					entry.geneId = g;
    					entry.transcriptIds.add(entry.transcript.tID);
    					entry.tool = InputConvertorConstants.REF_HEADER_ID;
    				}
    				fastaEntries.addAll(entries);
    			} else if(!isOnlyPC) {
    				ArrayList<FastaEntry> entries = this.getTranslation(t, isOnlyPC);
    				// put gene ID
    				for(FastaEntry entry : entries) {
    					entry.geneId = g;
    					entry.transcriptIds.add(entry.transcript.tID);
    					entry.tool = InputConvertorConstants.EXON_TRANSLATION_HEADER_ID;
    				}
    				fastaEntries.addAll(entries);
    			}
    			
        	}
        	
        });
        
        // remove duplications
        return FastaEntry.removeDuplications(fastaEntries);
	}
}
