package zhanglab.inputconvertor.data;

import java.util.ArrayList;
import java.util.Hashtable;

import zhanglab.inputconvertor.function.Translator;

public class FastaEntry {

	public String tool;
	public String sequence;
	public String header;
	public ArrayList<String> descriptions = new ArrayList<String>();
	
	// TODO
	public String table;
	
	public static ArrayList<FastaEntry> enumerateFastaEntry (Transcript t, ArrayList<Exon> exons) {
		ArrayList<FastaEntry> entries = new ArrayList<FastaEntry>();
		String strand = t.strand;
		Hashtable<String, String> indelHash = new Hashtable<String, String>();
		for(Exon exon : exons) {
			for(Mutation ins :exon.inss) {
				indelHash.put(ins.key, "");
			}
			for(Mutation del :exon.dels) {
				indelHash.put(del.key, "");
			}
		}
		
		StringBuilder wildSeq = new StringBuilder();
		for(Exon exon : exons) {
			wildSeq.append(exon.nucleotide);
		}
		
		// apply germline SNPs
		StringBuilder sequenceWithGermSNP = new StringBuilder();
		for(Exon exon : exons) {
			sequenceWithGermSNP.append(exon.getSequenceWithSNPs(true, false));
		}
		
		StringBuilder sequenceWithSomaSNP = new StringBuilder();
		for(Exon exon : exons) {
			sequenceWithSomaSNP.append(exon.getSequenceWithSNPs(true, true));
		}
		
		indelHash.forEach((key, nil)->{
			for(Exon exon : exons) {
				
			}
		});
		
		
		
		// wildtype translation
		for(int frame=0; frame<3; frame++) {
			if(strand.equalsIgnoreCase("+")) {
				Translator.translation(exons, frame);
			} else {
				Translator.reverseComplementTranslation(exons, frame);
			}
		}
		
		
		return entries;
	}
}
