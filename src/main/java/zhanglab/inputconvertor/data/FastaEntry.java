package zhanglab.inputconvertor.data;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;

import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.function.Translator;

public class FastaEntry {

	public int idx = -1;
	public String tool;
	public String sequence;
	public Transcript transcript;
	public String geneId;
	public ArrayList<String> transcriptIds = new ArrayList<String>();
	public int frame = -1;
	public String description = null;
	public String originHeader = null;
	public String mutationMark = null;
	
	public String toHeader(String uniqueId) {
		if(mutationMark != null) {
			return this.tool+this.mutationMark+this.idx+"_"+uniqueId;
		}
		
		return this.tool+this.idx+"_"+uniqueId;
	}
	
	public String toMeta (String uniqueId) {
		if(originHeader != null) {
			return toHeader(uniqueId)+"\t"+originHeader.replace("|", "\t");
		}
		
		String transcriptId = "";
		for(String id : transcriptIds) {
			if(transcriptId.length() != 0) {
				transcriptId += ",";
			}
			transcriptId += id;
		}
		return toHeader(uniqueId) +"\t" + this.geneId+"\t"+transcriptId+"\t"+this.frame+"\t"+transcript.strand+"\t"+this.description;
	}
	
	public String getKey () {
		return geneId+"_"+description+"_"+sequence;
	}
	
	
	public static ArrayList<Exon> getMutationNodes (Exon[] exonGraph) {
		ArrayList<Exon> mutExons = new ArrayList<Exon>();
		
		Exon exon = exonGraph[0];
		
		LinkedList<Exon> stack = new LinkedList<Exon>();
		stack.add(exon);
		
		Hashtable<Exon, Boolean> isVisited = new Hashtable<Exon, Boolean>();
		while(!stack.isEmpty()) {
			exon = stack.pollLast();
			
			if(isVisited.get(exon) != null) {
				continue;
			}
			isVisited.put(exon, true);
			
			if(exon.type != InputConvertorConstants.WILD) {
				mutExons.add(exon);
			}
			
			for(Exon nExon : exon.nextExons) {
				stack.add(nExon);
			}
		}
		
		return mutExons;
	}
	
	private static ArrayList<ArrayList<Exon>> rightEntries (Exon exon) {
		int MAX_FLANK_AA_SIZE = 14;
		int limit = MAX_FLANK_AA_SIZE * 3 + 2;
		
		ArrayList<ArrayList<Exon>> mutEntries = new ArrayList<ArrayList<Exon>>();
		
		LinkedList<Exon> exonStack = new LinkedList<Exon>();
		LinkedList<Exon> paths = new LinkedList<Exon>();
		int length = 0;
		exonStack.add(exon);
		
		while(!exonStack.isEmpty()) {
			Exon newExon = exonStack.pollLast();
			
			while(!paths.isEmpty()) {
				Exon path = paths.peekLast();
				boolean isRemoved = false;
				if(newExon.type == InputConvertorConstants.INS) {
					if(path.start > newExon.start) {
						paths.pollLast();
						isRemoved = true;
					}
				} else {
					if(path.start >= newExon.start) {
						paths.pollLast();
						isRemoved = true;
					}
				}
				
				if(isRemoved) {
					if(path.type == InputConvertorConstants.WILD) {
						length -= path.nucleotide.length();
					} else {
						length -= path.mutation.altSeq.length();
					}
				} else {
					break;
				}
			}
			
			if(newExon.start != Integer.MAX_VALUE) {
				if(newExon != exon) {
					paths.addLast(newExon);
				}
				
				if(newExon.type == InputConvertorConstants.WILD) {
					length += newExon.nucleotide.length();
				} else {
					length += newExon.mutation.altSeq.length();
				}
			}
			
			// reaching the end of transcript or exceeding length limitation
			if(length >= limit || newExon.start == Integer.MAX_VALUE) {
				ArrayList<Exon> mExons = new ArrayList<Exon>();
				Iterator<Exon> sExons = (Iterator<Exon>) paths.iterator();
				while(sExons.hasNext()) {
					mExons.add(sExons.next());
				}
				mutEntries.add(mExons);
			} else {
				for(Exon nExon : newExon.nextExons) {
					exonStack.add(nExon);
				}
			}
		}
		
		return mutEntries;
	}
	
	private static ArrayList<ArrayList<Exon>> leftEntries (Exon exon) {
		int MAX_FLANK_AA_SIZE = 14;
		int limit = MAX_FLANK_AA_SIZE * 3 + 2;
		
		ArrayList<ArrayList<Exon>> mutEntries = new ArrayList<ArrayList<Exon>>();
		
		LinkedList<Exon> exonStack = new LinkedList<Exon>();
		LinkedList<Exon> paths = new LinkedList<Exon>();
		int length = 0;
		exonStack.add(exon);
		
		while(!exonStack.isEmpty()) {
			Exon newExon = exonStack.pollLast();
			
			while(!paths.isEmpty()) {
				Exon path = paths.peekFirst();
				boolean isRemoved = false;
				if(newExon.type == InputConvertorConstants.INS) {
					if(path.end < newExon.end) {
						paths.pollFirst();
						isRemoved = true;
					}
				} else {
					if(path.end <= newExon.end) {
						paths.pollFirst();
						isRemoved = true;
					}
				}
				
				if(isRemoved) {
					if(path.type == InputConvertorConstants.WILD) {
						length -= path.nucleotide.length();
					} else {
						length -= path.mutation.altSeq.length();
					}
				} else {
					break;
				}
			}
			
			
			if(newExon.start != -1) {
				if(newExon != exon) {
					paths.addFirst(newExon);
				}
				
				if(newExon.type == InputConvertorConstants.WILD) {
					length += newExon.nucleotide.length();
				} else {
					length += newExon.mutation.altSeq.length();
				}
			}
			
			// reaching the start of transcript or exceeding length limitation
			if(length >= limit || newExon.start == -1) {
				ArrayList<Exon> mExons = new ArrayList<Exon>();
				Iterator<Exon> sExons = (Iterator<Exon>) paths.iterator();
				while(sExons.hasNext()) {
					mExons.add(sExons.next());
				}
				mutEntries.add(mExons);
			} else {
				for(Exon pExon : newExon.prevExons) {
					exonStack.add(pExon);
				}
			}
			
		}
		
		return mutEntries;
	}
	

	public static ArrayList<FastaEntry> enumerateCombinatorialMutations (Exon[] exonGraph) {
		ArrayList<FastaEntry> mutEntries = new ArrayList<FastaEntry>();
		// find mutation exons
		ArrayList<Exon> mutExons = getMutationNodes(exonGraph);
		
		for(Exon mutExon : mutExons) {
			ArrayList<ArrayList<Exon>> leftExons = leftEntries(mutExon);
			ArrayList<ArrayList<Exon>> rightExons = rightEntries(mutExon);
			
			int leftSize = leftExons.size();
			int rightSize = rightExons.size();
			
			for(int i=0; i<leftSize; i++) {
				ArrayList<Exon> leftExon = leftExons.get(i);
				String leftSeq = getSequenceFromExonList(leftExon);
				for(int j=0; j<rightSize; j++) {
					ArrayList<Exon> rightExon = rightExons.get(j);
					String rightSeq = getSequenceFromExonList(rightExon);
					
					FastaEntry entry = new FastaEntry();
					entry.description = getMutationDescription(leftExon) + getMutationDescription(mutExon) + getMutationDescription(rightExon);
					entry.sequence = leftSeq+mutExon.mutation.altSeq+rightSeq;
					entry.mutationMark = InputConvertorConstants.MUTATION_HEADER_ID;
					mutEntries.add(entry);
				}
			}
 		}
		
		return mutEntries;
	}
	
	public static String getSequenceFromExonList (ArrayList<Exon> exons) {
		StringBuilder sequence = new StringBuilder();
		
		for(Exon exon : exons) {
			if(exon.type == InputConvertorConstants.WILD) {
				sequence.append(exon.nucleotide);
			} else {
				sequence.append(exon.mutation.altSeq);
			}
		}
		
		return sequence.toString();
	}
	
	public static ArrayList<FastaEntry> enumerateFastaEntry (GenomeLoader refGenome, Transcript t, Collection<Exon> exons) {
		Exon[] exonGraph = Exon.toExonGraph(exons);
		refGenome.setSequence(t.chr, exonGraph);
		ArrayList<FastaEntry> entries = new ArrayList<FastaEntry>();
		
		ArrayList<FastaEntry> ntEntries = enumerateCombinatorialMutations(exonGraph);
		
		
		FastaEntry refEntry = new FastaEntry();
		
		StringBuilder seq = new StringBuilder();
		Exon exon = exonGraph[0];
		
		while(exon.start != Integer.MAX_VALUE) {
			if(exon.nucleotide != null) {
				seq.append(exon.nucleotide);
			}
			
			for(Exon nExon : exon.nextExons) {
				if(nExon.type == InputConvertorConstants.WILD) {
					exon = nExon;
					break;
				}
			}
		}
		refEntry.sequence = seq.toString();
		refEntry.description = getMutationDescription(exons);
		
		ntEntries.add(refEntry);
		
		for(FastaEntry fastaEntry : ntEntries) {
			fastaEntry.transcript = t;
			
			for(int frame=0; frame<3; frame++) {
				FastaEntry aaEntry = new FastaEntry();
				String peptide = null;
				if(t.strand.equalsIgnoreCase("+")) {
					peptide = Translator.translation(fastaEntry.sequence, frame);
				} else {
					peptide = Translator.reverseComplementTranslation(fastaEntry.sequence, frame);
				}
				aaEntry.description = fastaEntry.description;
				aaEntry.frame = frame;
				aaEntry.sequence = peptide;
				aaEntry.transcript = t;
				aaEntry.mutationMark = fastaEntry.mutationMark;
				
				entries.add(aaEntry);
			}
			
		}
		
		return entries;
	}
	
	public static ArrayList<FastaEntry> enumerateFastaEntryCDS (GenomeLoader refGenome, Transcript t, Collection<Exon> exons) {
		Exon[] exonGraph = Exon.toExonGraph(exons);
		refGenome.setSequence(t.chr, exonGraph);
		ArrayList<FastaEntry> entries = new ArrayList<FastaEntry>();
		
		FastaEntry refEntry = new FastaEntry();
		
		StringBuilder seq = new StringBuilder();
		Exon exon = exonGraph[0];
		
		while(exon.start != Integer.MAX_VALUE) {
			if(exon.nucleotide != null) {
				seq.append(exon.nucleotide);
			}
			
			for(Exon nExon : exon.nextExons) {
				if(nExon.type == InputConvertorConstants.WILD) {
					exon = nExon;
					break;
				}
			}
		}
		refEntry.sequence = seq.toString();
		refEntry.description = getMutationDescription(exons);
		refEntry.transcript = t;
		
		int longestFrameIdx = 0;
		int longestLength = 0;
		
		// determine longest frame (to deal with some fragile RNAs)
		for(int frame=0; frame<3; frame++) {
			String peptide = null;
			if(t.strand.equalsIgnoreCase("+")) {
				peptide = Translator.translation(refEntry.sequence, frame);
			} else {
				peptide = Translator.reverseComplementTranslation(refEntry.sequence, frame);
			}
			
			String[] xPeptides = peptide.split("X");
			for(int i=0; i<xPeptides.length; i++) {
				if(longestLength < xPeptides[i].length()) {
					longestLength = xPeptides[i].length();
					longestFrameIdx = frame;
				}
			}
		}
		
		FastaEntry aaEntry = new FastaEntry();
		String peptide = null;
		if(t.strand.equalsIgnoreCase("+")) {
			peptide = Translator.translation(refEntry.sequence, longestFrameIdx);
		} else {
			peptide = Translator.reverseComplementTranslation(refEntry.sequence, longestFrameIdx);
		}
		aaEntry.description = refEntry.description;
		aaEntry.frame = longestFrameIdx;
		aaEntry.sequence = peptide;
		aaEntry.transcript = t;
		aaEntry.mutationMark = refEntry.mutationMark; // must be null
		
		entries.add(aaEntry);
		
		
		return entries;
	}

	private static String getMutationDescription (Exon exon) {
		StringBuilder desc = new StringBuilder();
		
		desc.append("@").append(exon.getMutationDescription());
		
		return desc.toString();
	}
	
	private static String getMutationDescription (Collection<Exon> exons) {
		StringBuilder desc = new StringBuilder();
		
		for(Exon exon : exons) {
			desc.append("@").append(exon.getMutationDescription());
		}
		
		return desc.toString();
	}
	
	public static ArrayList<FastaEntry> removeDuplications (ArrayList<FastaEntry> entries) {
 		// it is possible to appear duplicated peptides because of several reference transcripts. 
 		Hashtable<String, FastaEntry> removeDuplication = new Hashtable<String, FastaEntry>();
 		ArrayList<FastaEntry> uniqueFastaEntries = new ArrayList<FastaEntry>();
 		int idx = 1;
         for(FastaEntry entry : entries) {
        	 FastaEntry firstEntry = removeDuplication.get(entry.getKey());
         	if(firstEntry != null) {
         		for(String tID1 : entry.transcriptIds) {
         			boolean isIncluded = false;
         			for(String tID2 : firstEntry.transcriptIds) {
         				if(tID1.equalsIgnoreCase(tID2)) {
         					isIncluded = true;
         				}
         			}
         			if(!isIncluded) {
         				firstEntry.transcriptIds.add(tID1);
         			}
         		}
         		continue;
         	} else {
         		removeDuplication.put(entry.getKey(), entry);
             	entry.idx = idx++;
             	uniqueFastaEntries.add(entry);
         	}
         }
         
         return uniqueFastaEntries;
	}
}
