package zhanglab.inputconvertor.input;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
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
	public String logFileName = null;
	
	public StringTie (File file) throws IOException {
		System.out.println("## StringTie ##");
		System.out.println("Load "+file.getName());
		this.gtf = new GTFLoader(file);
		
		this.logFileName = file.getAbsolutePath()+".SR.log";
		System.out.println("## StringTie ##");
		System.out.println("Generate log file... "+this.logFileName);
		this.writeLog();
	}
	
	public void enrollGenomeSequence (GenomeLoader refGenome) {
		System.out.println("## StringTie ##");
		System.out.println("Enroll reference genome");
		this.refGenome = refGenome;
	}
	
	private void writeLog () throws IOException {
		BufferedWriter BW = new BufferedWriter(new FileWriter(this.logFileName));
		
		int[] nums = new int[11]; 
		nums[0] = this.gtf.geneToTranscripts.size();
		// nums[1]: total number of transcripts
		// nums[2]: a number of transcripts log2FPKM [0, 1)
		// nums[3]: a number of transcripts log2FPKM [1, 2)
		// nums[4]: a number of transcripts log2FPKM [2, 3)
		// nums[5]: a number of transcripts log2FPKM [3, 4)
		// nums[6]: a number of transcripts log2FPKM [4, 5)
		// nums[7]: a number of transcripts log2FPKM [5, 6)
		// nums[8]: a number of transcripts log2FPKM [6, 7)
		// nums[9]: a number of transcripts log2FPKM [7, 8)
		// nums[10]: a number of transcripts log2FPKM [8, 9)
		// nums[11]: a number of transcripts log2FPKM [9, 10)
		// nums[12]: a number of transcripts log2FPKM [10, 11)
		// nums[13]: a number of transcripts log2FPKM [11, 12)
		// nums[14]: a number of transcripts log2FPKM [12, 13)
		// nums[15]: a number of transcripts log2FPKM [13, 14)
		// nums[16]: a number of transcripts log2FPKM [14, ~]
		this.gtf.geneToTranscripts.forEach((g, ts)->{
			nums[1] += ts.size();
			for(Transcript t : ts) {
				double fpkm = t.FPKM;
				double log2FPKM = Math.log(fpkm+1)/Math.log(2);
				int idx = (int) log2FPKM;
				idx = idx + 2;
				
				if(idx >= nums.length) {
					idx = nums.length-1;
				}
				nums[idx]++;
			}
		});
		
		
		BW.append("Number of genes\t"+nums[0]);
		BW.newLine();
		BW.append("Number of transcripts\t"+nums[1]);
		BW.newLine();
		BW.append("Log2(FPKM+1)");
		BW.newLine();
		for(int i=2; i<nums.length; i++) {
			int min = i-2;
			int max = i-1;
			if(i == nums.length-1) {
				BW.append("["+min+",inf)\t"+nums[i]);
			} else {
				BW.append("["+min+","+max+")\t"+nums[i]);
			}
			
			BW.newLine();
		}
		
		BW.close();
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
    			if(t.FPKM < FPKMthreshold) {
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
