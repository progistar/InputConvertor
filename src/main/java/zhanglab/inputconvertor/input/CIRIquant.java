package zhanglab.inputconvertor.input;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import zhanglab.inputconvertor.data.Exon;
import zhanglab.inputconvertor.data.FastaEntry;
import zhanglab.inputconvertor.data.GTFLoader;
import zhanglab.inputconvertor.data.GenomeLoader;
import zhanglab.inputconvertor.data.Transcript;
import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.function.Translator;

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
	
	private String getTranslation (ArrayList<Exon> exons, String chr, int tStart, int tEnd, int frame, String strand) {
		StringBuilder protein = new StringBuilder();
		
		ArrayList<Exon> nExons = new ArrayList<Exon>();
		boolean hasStart = false;
		boolean hasEnd = false;
		for(Exon e : exons) {
			if(e.start <= tStart && e.end >= tStart) {
				hasStart = true;
			}
			
			if(hasStart && !hasEnd) {
				nExons.add(new Exon(e.start, e.end));
			}
			
			if(e.start <= tEnd && e.end >= tEnd) {
				hasEnd = true;
			}
		}
		
		if(hasStart && hasEnd) {
    		for(int i=0; i<nExons.size(); i++) {
    			Exon exon = nExons.get(i);
    			if(i == 0) {
    				exon.start = tStart;
    			} 
    			if(i == nExons.size()-1) {
    				exon.end = tEnd;
    			}
    		}
    		refGenome.setSequence(chr, nExons);
    		int len = getLengthOfExons(nExons);
    		
    		// 2 x linear sequences
    		nExons.addAll(nExons);
    		
			// if the sequence is less than 30 nt, then we do more concatenate the sequence.
			if(len < MAX_FLANK_AA_SIZE * 3) {
				nExons.addAll(nExons);
			}
			
			if(strand.equalsIgnoreCase("+")) {
    			protein.append(Translator.translation(nExons, frame));
    		} else {
    			protein.append(Translator.reverseComplementTranslation(nExons, frame));
    		}
			
			// max +-50AA
			int mid = protein.length() / 2;
			int leftMax = mid - MAX_FLANK_AA_SIZE;
			int rightMax = mid + MAX_FLANK_AA_SIZE;
			
			if(leftMax < 0) {
				leftMax = 0;
			}
			if(rightMax > protein.length()) {
				rightMax = protein.length();
			}
			
			StringBuilder leftSeq = new StringBuilder();
			for(int j=mid-1; j>=leftMax; j--) {
				if(protein.charAt(j) == 'X') {
					break;
				} else {
					leftSeq.append(protein.charAt(j));
				}
			}
			leftSeq = leftSeq.reverse();
			
			StringBuilder rightSeq = new StringBuilder();
			for(int j=mid; j<rightMax; j++) {
				if(protein.charAt(j) == 'X') {
					break;
				} else {
					rightSeq.append(protein.charAt(j));
				}
			}
			
			protein.setLength(0);
			
			if(leftSeq.length() > 0 && rightSeq.length() > 0) {
				protein.append(leftSeq.toString()).append(rightSeq.toString());
			}
		}
		
		return protein.toString();
	}
	
	public ArrayList<FastaEntry> getFastaEntry () throws IOException {
		ArrayList<FastaEntry> fastaEntries = new ArrayList<FastaEntry>();
		
        this.gtf.geneToTranscripts.forEach((g, ts)->{
    		ArrayList<Transcript> refTs = this.refGTF.geneToTranscripts.get(g);
    		
    		for(Transcript t : ts) {
    			int tStart = Integer.parseInt(t.start);
    			int tEnd = Integer.parseInt(t.end);
    			// it is possible to appear duplicated peptides because of several reference transcripts. 
    			Hashtable<String, String> removeDuplication = new Hashtable<String, String>();
    			
    			ArrayList<ArrayList<Exon>> exons = new ArrayList<ArrayList<Exon>>();
    			String protein = null;
    			
    			// intergenic
    			if(refTs == null) {
    				ArrayList<Exon> nExons = new ArrayList<Exon>();
    				nExons.add(new Exon(tStart, tEnd));
    				
    			} else {
        			// start translating
        			for(Transcript refT : refTs) {
        				exons.add(refT.exons);
        			}
        			// end 
    			}
    			
    			for(ArrayList<Exon> nExons : exons) {
    				for(int frame = 0; frame < 3; frame++) {
    					protein = this.getTranslation(nExons, t.chr, tStart, tEnd, frame, t.strand);
    					
    					String[] proteins = protein.split("X");
        				int idx = 1;
        				for(int i=0; i<proteins.length; i++) {
        					protein = proteins[i];
        					
        					// cut
        					if(protein.length() < InputConvertorConstants.MIN_PEPT_LEN) continue;
        					
        					if(removeDuplication.get(protein) != null) continue;
        					removeDuplication.put(protein, "");
        					
    						FastaEntry entry = new FastaEntry();
    						entry.tool = InputConvertorConstants.CIRIQUANT_HEADER_ID;
    						entry.header = t.tID+"|"+frame+"|"+idx++;
    						entry.sequence = protein;
    						fastaEntries.add(entry);
        				}
    				}
    			}
        	}
        	
        });
        
        return fastaEntries;
	}
	
	public int getLengthOfExons (ArrayList<Exon> exons) {
		int len = 0;
		for(Exon exon : exons ) {
			len += exon.nucleotide.length();
		}
		
		return len;
	}
}
