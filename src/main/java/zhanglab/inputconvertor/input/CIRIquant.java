package zhanglab.inputconvertor.input;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import zhanglab.inputconvertor.data.Exon;
import zhanglab.inputconvertor.data.GTFLoader;
import zhanglab.inputconvertor.data.GenomeLoader;
import zhanglab.inputconvertor.data.Transcript;
import zhanglab.inputconvertor.function.Translator;

public class CIRIquant {

	public static final int MAX_FLANK_SIZE = 14;
	
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
		StringBuilder nucleotide = new StringBuilder();
		StringBuilder protein = new StringBuilder();
		
		ArrayList<Exon> nExons = new ArrayList<Exon>();
		boolean hasStart = false;
		boolean hasEnd = false;
		for(Exon e : exons) {
			if(e.start <= tStart && e.end >= tStart) {
				hasStart = true;
			}
			
			if(hasStart && !hasEnd) {
				nExons.add(e);
			}
			
			if(e.start <= tEnd && e.end >= tEnd) {
				hasEnd = true;
			}
		}
		
		if(hasStart && hasEnd) {
    		for(int i=0; i<nExons.size(); i++) {
    			int eStart = nExons.get(i).start - 1;
    			int eEnd = nExons.get(i).end;
    			if(i == 0) {
    				eStart = tStart - 1;
    				
    			} 
    			if(i == nExons.size()-1) {
    				eEnd = tEnd;
    			}
    			nucleotide.append(refGenome.getSequence(chr, eStart, eEnd));
    		}
    		
    		// 2 x linear sequences
    		String originNucleotide = nucleotide.toString();
    		if(originNucleotide.length() <= 30) {
    			System.out.println(chr+":"+tStart+"|"+tEnd);
    			System.out.println(originNucleotide);
    		}
			nucleotide.append(nucleotide.toString());
			// if the sequence is less than 300 nt, then we do more concatenate the sequence.
			if(nucleotide.length() < MAX_FLANK_SIZE * 3) {
				nucleotide.append(originNucleotide).append(originNucleotide);
			}
			
			if(strand.equalsIgnoreCase("+")) {
    			protein.append(Translator.translation(nucleotide.toString(), frame));
    		} else {
    			protein.append(Translator.reverseComplementTranslation(nucleotide.toString(), frame));
    		}
			
			// max +-50AA
			int mid = protein.length() / 2;
			int leftMax = mid - MAX_FLANK_SIZE;
			int rightMax = mid + MAX_FLANK_SIZE;
			
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
	
	public void writeEntry (BufferedWriter mainFasta, BufferedWriter mapper) throws IOException {
		
		
		
        this.gtf.geneToTranscripts.forEach((g, ts)->{
    		ArrayList<Transcript> refTs = this.refGTF.geneToTranscripts.get(g);
    		
    		for(Transcript t : ts) {
    			int tStart = Integer.parseInt(t.start);
    			int tEnd = Integer.parseInt(t.end);
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
    					if(protein.length() != 0) {
    						try {
    							mainFasta.append(">").append(t.tID).append("|").append(""+frame);
    							mainFasta.newLine();
    							mainFasta.append(protein.toString());
    							mainFasta.newLine();
    							
    						}catch(IOException ioe) {
    							
    						}
    					}
    				}
    			}
        	}
        	
        });
	}
}
