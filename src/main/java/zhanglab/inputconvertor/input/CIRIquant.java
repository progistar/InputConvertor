package zhanglab.inputconvertor.input;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;

import zhanglab.inputconvertor.data.Exon;
import zhanglab.inputconvertor.data.FastaEntry;
import zhanglab.inputconvertor.data.GTFLoader;
import zhanglab.inputconvertor.data.GenomeLoader;
import zhanglab.inputconvertor.data.Transcript;
import zhanglab.inputconvertor.env.InputConvertorConstants;

public class CIRIquant {

	public static final int MAX_FLANK_AA_SIZE = 14;
	
	public GTFLoader gtf = null;
	public GTFLoader refGTF = null;
	public GenomeLoader refGenome = null;
	
	public CIRIquant (File file) throws IOException {
		System.out.println("## CIRIquant ##");
		System.out.println("Load "+file.getName());
		this.gtf = new GTFLoader(file);
	}
	
	public void enrollReferenceGTF (GTFLoader refGTF) {
		System.out.println("## CIRIquant ##");
		System.out.println("Enroll reference GTF");
		this.refGTF = refGTF;
	}
	
	public void enrollGenomeSequence (GenomeLoader refGenome) {
		System.out.println("## CIRIquant ##");
		System.out.println("Enroll reference genome");
		this.refGenome = refGenome;
	}
	
	private ArrayList<FastaEntry> getTranslation (Transcript t, ArrayList<Exon> exons, int cStart, int cEnd) {
		ArrayList<FastaEntry> entries = new ArrayList<FastaEntry>();
		ArrayList<Exon> nExons = new ArrayList<Exon>();
		for(Exon e : exons) {
			if(e.end >= cStart && e.start <= cEnd) {
				Exon nExon = new Exon(t.chr, e.start, e.end);
				nExons.add(nExon);
			}
		}
		
		if(nExons.size() != 0) {
			int limit = MAX_FLANK_AA_SIZE * 3 + 2;
			LinkedList<Exon> circExons = new LinkedList<Exon>();
			// adjust start and end
			// :: intron-compatibility
			nExons.get(0).start = cStart;
			nExons.get(nExons.size()-1).end = cEnd;
			
			int leftExonIdx = nExons.size()-1;
			int rightExonIdx = leftExonIdx+1;
			
			nExons.addAll(nExons);
			
			while(getLengthOfExons(nExons) < 2 * (limit)) {
				leftExonIdx = nExons.size()-1;
				rightExonIdx = leftExonIdx+1;
				nExons.addAll(nExons);
			}
			
			// left shrink
			int leftSize = 0;
			for(int i=leftExonIdx; i>=0; i--) {
				Exon exon = nExons.get(i);
				Exon nExon = new Exon(t.chr, exon.start, exon.end);
				circExons.addFirst(nExon);
    			if(nExon.end - nExon.start + 1 + leftSize >= limit) {
    				nExon.start = nExon.end + 1 + leftSize - limit;
    				leftSize += (nExon.end - nExon.start + 1);
    				break;
    			}
    			leftSize += (nExon.end - nExon.start + 1);
			}
			
			// right shrink
			int rightSize = 0;
			for(int i=rightExonIdx; i<nExons.size(); i++) {
				Exon exon = nExons.get(i);
				Exon nExon = new Exon(t.chr, exon.start, exon.end);
				circExons.add(nExon);
    			if(nExon.end - nExon.start + 1 + rightSize >= limit) {
    				nExon.end = limit + nExon.start -1 - rightSize;
    				rightSize += (nExon.end - nExon.start + 1);
    				break;
    			}
    			rightSize += (nExon.end - nExon.start + 1);
			}
			
			//
			entries = FastaEntry.enumerateFastaEntry(refGenome, t, circExons, 1);
		}
		return entries;
	}
	
	public ArrayList<FastaEntry> getFastaEntry () throws IOException {
		ArrayList<FastaEntry> fastaEntries = new ArrayList<FastaEntry>();
		
        this.gtf.geneToTranscripts.forEach((g, ts)->{
    		ArrayList<Transcript> refTs = this.refGTF.geneToTranscripts.get(g);
    		
    		
    		for(Transcript t : ts) {
    			int tStart = Integer.parseInt(t.start);
    			int tEnd = Integer.parseInt(t.end);
    			// it is possible to appear duplicated peptides because of several reference transcripts. 
    			ArrayList<ArrayList<Exon>> exons = new ArrayList<ArrayList<Exon>>();
    			ArrayList<Transcript> matchedTranscripts = new ArrayList<Transcript>();
    			
    			
    			// intergenic
    			if(refTs == null) {
    				ArrayList<Exon> nExons = new ArrayList<Exon>();
    				nExons.add(new Exon(t.chr, tStart, tEnd));
    				matchedTranscripts.add(t);
    			} else {
        			// start translating
        			for(Transcript refT : refTs) {
        				exons.add(refT.exons);
        				matchedTranscripts.add(refT);
        			}
        			// end 
    			}
    			
    			// reflect reference transcript ID
    			for(int i=0; i<exons.size(); i++) {
    				ArrayList<Exon> nExons = exons.get(i);
    				Transcript matchedTranscript = matchedTranscripts.get(i);
    				ArrayList<FastaEntry> entries = this.getTranslation(matchedTranscript, nExons, tStart, tEnd);
    				// put gene ID
    				for(FastaEntry entry : entries) {
    					String geneName = this.refGTF.geneToGeneName.get(g);
    					entry.geneName = geneName == null ? g : geneName;
    					entry.geneId = g;
    					entry.transcriptIds.add(entry.transcript.tID);
    					entry.tool = InputConvertorConstants.CIRIQUANT_HEADER_ID;
    				}
    				fastaEntries.addAll(entries);
    			}
        	}
        	
        });
        
        // remove duplications
        return FastaEntry.removeDuplications(fastaEntries);
	}
	
	public int getLengthOfExons (ArrayList<Exon> exons) {
		int len = 0;
		for(Exon exon : exons ) {
			len += (exon.end - exon.start + 1);
		}
		
		return len;
	}
}
