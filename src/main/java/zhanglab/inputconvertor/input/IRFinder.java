package zhanglab.inputconvertor.input;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.LinkedList;

import zhanglab.inputconvertor.data.Exon;
import zhanglab.inputconvertor.data.FastaEntry;
import zhanglab.inputconvertor.data.GTFLoader;
import zhanglab.inputconvertor.data.GenomeLoader;
import zhanglab.inputconvertor.data.Transcript;
import zhanglab.inputconvertor.env.InputConvertorConstants;

public class IRFinder {
	public static final int MAX_FLANK_AA_SIZE = 14;
	public Hashtable<String, ArrayList<Transcript>> geneToTranscripts = new Hashtable<String, ArrayList<Transcript>>();
	public File file = null;
	public GTFLoader refGTF = null;
	public GenomeLoader refGenome = null;
	
	public void clear() {
		geneToTranscripts.clear();
	}
	
	public IRFinder (File file) throws IOException {
		System.out.println("## IRFinder ##");
		System.out.println("Load "+file.getName());
		this.file = file;
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		String[] header = BR.readLine().split("\t");
		
		int chrIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.IRFINDER_CHR_FIELD_NAME);
		int startIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.IRFINDER_START_FIELD_NAME);
		int endIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.IRFINDER_END_FIELD_NAME);
		int nameIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.IRFINDER_NAME_FIELD_NAME);
		int strandIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.IRFINDER_STRAND_FIELD_NAME);
		int warningIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.IRFINDER_WARNINGS_FIELD_NAME);
		
		Hashtable<String, Integer> summary = new Hashtable<String, Integer>();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String name = fields[nameIdx];
			String chr = fields[chrIdx];
			// IRFinder reports [start, end) with zero-based index
			int start = Integer.parseInt(fields[startIdx])+1;
			int end = Integer.parseInt(fields[endIdx]);
			//////////////////////////////////////////////////////
			String strand = fields[strandIdx];
			String warnings = fields[warningIdx];
			String geneId = name.split("\\/")[1];
			String id = chr+":"+start+"-"+end+"|"+strand+"|"+name;
			
			Integer cnt = summary.get(warnings);
			if(cnt == null) {
				cnt = 0;
			}
			summary.put(warnings, cnt+1);
			// pass warnings
			if(!warnings.equalsIgnoreCase("-")) {
				continue;
			}
			
			Transcript transcript = new Transcript();
			transcript.tID = id;
			transcript.start = String.valueOf(start);
			transcript.end = String.valueOf(end);
			transcript.chr = chr;
			transcript.strand = strand;
			transcript.attrs = warnings;
			Exon exon = new Exon(transcript.chr, start, end);
			transcript.exons.add(exon);
			
			ArrayList<Transcript> transcripts = geneToTranscripts.get(geneId);
			if(transcripts == null) {
				transcripts = new ArrayList<Transcript>();
				geneToTranscripts.put(geneId, transcripts);
			}
			
			transcripts.add(transcript);
		}
		
		BR.close();
		
		// summary
		summary.forEach((w, n)->{
			System.out.println(w+"\t"+n);
		});
		
	}
	
	public void enrollReferenceGTF (GTFLoader refGTF) {
		System.out.println("## IRFinder ##");
		System.out.println("Enroll reference GTF");
		this.refGTF = refGTF;
	}
	
	public void enrollGenomeSequence (GenomeLoader refGenome) {
		System.out.println("## IRFinder ##");
		System.out.println("Enroll reference genome");
		this.refGenome = refGenome;
	}
	
	private ArrayList<FastaEntry> getTranslation (Transcript t, ArrayList<Exon> exons, int iStart, int iEnd) {
		ArrayList<FastaEntry> entries = new ArrayList<FastaEntry>();
		int startExon = -1;
		int endExon = -1;
		for(int i=0; i<exons.size(); i++) {
			Exon e = exons.get(i);
			if(e.start <= iStart-1 && e.end >= iStart-1) {
				startExon = i;
			}
			if(e.start <= iEnd+1 && e.end >= iEnd+1) {
				endExon = i;
			}
		}
		
		Exon intron = new Exon(t.chr, iStart, iEnd);
		
		if(startExon != -1 && endExon != -1) {
			int limit = MAX_FLANK_AA_SIZE * 3 + 2;
			LinkedList<Exon> nExons = new LinkedList<Exon>();
			// left exons
			int leftSize = 0;
			for(int i=startExon; i>=0; i--) {
				Exon exon = exons.get(i);
				Exon nExon = new Exon(t.chr, exon.start, exon.end);
				nExons.addFirst(nExon);
    			if(nExon.end - nExon.start + 1 + leftSize >= limit) {
    				nExon.start = nExon.end + 1 + leftSize - limit;
    				leftSize += (nExon.end - nExon.start + 1);
    				break;
    			}
    			leftSize += (nExon.end - nExon.start + 1);
			}
			nExons.add(intron);
			// right exons
			int rightSize = 0;
			for(int i=endExon; i<exons.size(); i++) {
				Exon exon = exons.get(i);
				Exon nExon = new Exon(t.chr, exon.start, exon.end);
				nExons.add(nExon);
    			if(nExon.end - nExon.start + 1 + rightSize >= limit) {
    				nExon.end = limit + nExon.start -1 - rightSize;
    				rightSize += (nExon.end - nExon.start + 1);
    				break;
    			}
    			rightSize += (nExon.end - nExon.start + 1);
			}
			
			refGenome.setSequence(t.chr, nExons);
			entries = FastaEntry.enumerateFastaEntry(t, nExons);
		}
		return entries;
	}
	
	public ArrayList<FastaEntry> getFastaEntry () throws IOException {
		ArrayList<FastaEntry> fastaEntries = new ArrayList<FastaEntry>();
		
        this.geneToTranscripts.forEach((g, ts)->{
    		ArrayList<Transcript> refTs = this.refGTF.geneToTranscripts.get(g);
    		for(Transcript t : ts) {
    			int tStart = Integer.parseInt(t.start);
    			int tEnd = Integer.parseInt(t.end);
    			
    			ArrayList<ArrayList<Exon>> exons = new ArrayList<ArrayList<Exon>>();
    			// If there is no matched transcripts, it means that the intron is located in "intergenic" region.
    			// But! it is impossible.
    			// There is only case: the version is not matched.
    			if(refTs == null) {
    				System.out.println("Severe:: The reference transcript version is not matched to that of IRFinder");
    				
    			} else {
        			// start translating
        			for(Transcript refT : refTs) {
        				exons.add(refT.exons);
        			}
        			// end 
    			}
    			
    			for(ArrayList<Exon> nExons : exons) {
    				ArrayList<FastaEntry> entries = this.getTranslation(t, nExons, tStart, tEnd);
    				// put gene ID
    				for(FastaEntry entry : entries) {
    					entry.geneId = g;
    				}
    				fastaEntries.addAll(entries);
    			}
        	}
        	
        });
        
        // add tool info
		// it is possible to appear duplicated peptides because of several reference transcripts. 
		Hashtable<String, String> removeDuplication = new Hashtable<String, String>();
		ArrayList<FastaEntry> uniqueFastaEntries = new ArrayList<FastaEntry>();
		int idx = 1;
        for(FastaEntry entry : fastaEntries) {
        	if(removeDuplication.get(entry.getKey()) != null) {
        		continue;
        	}
        	removeDuplication.put(entry.getKey(), "");
        	entry.tool = InputConvertorConstants.IRFINDER_HEADER_ID;
        	entry.idx = idx++;
        	uniqueFastaEntries.add(entry);
        }
        
        fastaEntries.clear();
        
        return uniqueFastaEntries;
	}
}
