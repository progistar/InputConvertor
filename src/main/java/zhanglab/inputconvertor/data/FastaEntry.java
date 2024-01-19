package zhanglab.inputconvertor.data;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Hashtable;

import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.function.Translator;

public class FastaEntry {

	public int idx = -1;
	public String tool;
	public String sequence;
	public Transcript transcript;
	public String geneId;
	public int frame = -1;
	public String description = null;
	
	public String toHeader() {
		return  this.tool+this.idx+"|"+this.geneId+"|f:"+this.frame+"|s:"+transcript.strand+"|"+this.description;
	}
	
	public String getKey () {
		return geneId+"_"+description+"_"+sequence;
	}
	
	// TODO
	public String table;
	
	public static ArrayList<FastaEntry> enumerateFastaEntry (Transcript t, Collection<Exon> exons) {
		ArrayList<FastaEntry> entries = new ArrayList<FastaEntry>();
		//######### DEBUG FOR MUTANT ONLY ########################################################################################################
		boolean isRun = false;
		for(Exon exon : exons) {
			if(exon.isMutant()) {
				isRun = true;
			}
		}
		//if(!isRun) return entries;
		//######### DEBUG FOR MUTANT ONLY ########################################################################################################		
		
		Hashtable<String, Byte> indelHash = new Hashtable<String, Byte>();
		indelHash.put("wildtype", InputConvertorConstants.NON_INDEL);
		for(Exon exon : exons) {
			for(Mutation ins :exon.inss) {
				indelHash.put(ins.key, InputConvertorConstants.INS);
			}
			for(Mutation del :exon.dels) {
				indelHash.put(del.key, InputConvertorConstants.DEL);
			}
		}
		
		indelHash.forEach((indelID, type)->{
			StringBuilder wildSeqWithINDEL = new StringBuilder();
			StringBuilder germWithINDEL = new StringBuilder();
			StringBuilder germsomaWithINDEL = new StringBuilder();
			
			for(Exon exon : exons) {
				if(type == InputConvertorConstants.NON_INDEL || type == InputConvertorConstants.INS) {
					wildSeqWithINDEL.append(exon.getSequenceWithINS(indelID, false, false));
					germWithINDEL.append(exon.getSequenceWithINS(indelID, true, false));
					germsomaWithINDEL.append(exon.getSequenceWithINS(indelID, true, true));
				} else if(type == InputConvertorConstants.NON_INDEL ||type == InputConvertorConstants.DEL) {
					wildSeqWithINDEL.append(exon.getSequenceWithDEL(indelID, false, false));
					germWithINDEL.append(exon.getSequenceWithDEL(indelID, true, false));
					germsomaWithINDEL.append(exon.getSequenceWithDEL(indelID, true, true));
				}
			}
			
			FastaEntry entry = new FastaEntry();
			entry.description = getMutationDescription(exons, indelID, false, false);
			entry.sequence = wildSeqWithINDEL.toString();
			entry.transcript = t;
			entries.add(entry);
			
			entry = new FastaEntry();
			entry.description = getMutationDescription(exons, indelID, true, false);
			entry.sequence = germWithINDEL.toString();
			entry.transcript = t;
			entries.add(entry);
			
			entry = new FastaEntry();
			entry.description = getMutationDescription(exons, indelID, true, true);
			entry.sequence = germsomaWithINDEL.toString();
			entry.transcript = t;
			entries.add(entry);
			
		});
		
		// remove duplicated sequences by exon information
		Hashtable<String, FastaEntry> uniqueEntries = new Hashtable<String, FastaEntry>();
		for(FastaEntry entry : entries) {
			uniqueEntries.put(entry.description, entry);
		}
		entries.clear();
		
		// translation
		uniqueEntries.forEach((desc, entry)->{
			for(int frame=0; frame<3; frame++) {
				FastaEntry aaEntry = new FastaEntry();
				String peptide = null;
				if(t.strand.equalsIgnoreCase("+")) {
					peptide = Translator.translation(entry.sequence, frame);
				} else {
					peptide = Translator.reverseComplementTranslation(entry.sequence, frame);
				}
				
				aaEntry.description = entry.description;
				aaEntry.frame = frame;
				aaEntry.sequence = peptide;
				aaEntry.transcript = t;
				
				entries.add(aaEntry);
			}
		});
		
		return entries;
	}

	private static String getMutationDescription (Collection<Exon> exons, String indelID, boolean germSNP, boolean somaSNP) {
		StringBuilder desc = new StringBuilder();
		
		for(Exon exon : exons) {
			desc.append("@").append(exon.getMutationDescription(indelID, germSNP, somaSNP));
		}
		
		return desc.toString();
	}
}
