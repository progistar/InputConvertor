package zhanglab.inputconvertor.data;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.function.Translator;
import zhanglab.inputconvertor.run.MakeNeoDB;

public class FastaEntry {

	public int idx = -1;
	public String tool;
	public String sequence;
	public String nucleotide;
	public Transcript transcript;
	public String transcriptId;
	public String geneId;
	public String geneName;
	public String strand;
	public String frame = ".";
	public String description = null;
	public String originHeader = null;
	public String mutationMark = null;
	
	public String toHeader(String uniqueId) {
		if(uniqueId == null) {
			return this.tool;
		}
		return this.tool+"_"+uniqueId;
	}
	
	public String toMeta (String uniqueId) {
		if(originHeader != null) {
			return toHeader(uniqueId)+"|"+originHeader;
		}
		while(this.description.startsWith("@")) {
			this.description = this.description.substring(1);
		}
		
		if(this.mutationMark == null) {
			this.mutationMark = "wild";
		}
		
		return toHeader(uniqueId) +"|" + 
		this.transcriptId + "#" +
		this.idx+ 
		"|"+this.geneId+" GN="+this.geneName+" FR="+this.frame+" SR="+this.strand+" gene_site="+this.description.replace("@"," ");
	}
	
	public String toMetaSimple () {
		while(this.description.startsWith("@")) {
			this.description = this.description.substring(1);
		}
		
		if(this.mutationMark == null) {
			this.mutationMark = "wild";
		}
		
		String[] locations = this.description.split("@");
		String startSite = locations[0].split("\\-")[0];
		// put the strand info
		startSite = startSite.replace(":", transcript.strand+":");
		String endSite = locations[locations.length-1].split("\\-")[1];
		
		
		StringBuilder mutations = new StringBuilder();
		for(String location : locations) {
			if(location.contains("[SNP]") ||
					location.contains("[MNP]") ||
					location.contains("[INS]") ||
					location.contains("[DEL]")) {
				if(mutations.length() > 0) {
					mutations.append(" ");
				}
				mutations.append(location);
			}
		}
		
		return this.tool +"|" + 
		transcriptId + "_" + 
		startSite+ "-" +endSite + "_" +
		this.frame + "_"+ 
		this.mutationMark+"#"+ this.idx+ 
		"|"+this.geneId+" "+this.transcript.classCode+" "+mutations.toString();
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
		int limit = MakeNeoDB.maxFlankLength * 3;
		
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
				
				if(length >= limit) {
					paths.pollLast();
					isRemoved = true;
				} else if(newExon.type == InputConvertorConstants.INS) {
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
						length -= path.refNucleotide.length();
					} else {
						length -= path.altNucleotide.length();
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
					length += newExon.refNucleotide.length();
				} else {
					length += newExon.altNucleotide.length();
				}
			}
			
			// reaching the end of transcript or exceeding length limitation
			if(length >= limit || newExon.start == Integer.MAX_VALUE) {
				ArrayList<Exon> mExons = new ArrayList<Exon>();
				Iterator<Exon> sExons = (Iterator<Exon>) paths.iterator();
				int size = paths.size();
				while(sExons.hasNext()) {
					Exon e = sExons.next();
					if(--size == 0 && newExon.start != Integer.MAX_VALUE) {
						e = e.copyExon(e.start, e.end - (length - limit));
					}
					mExons.add(e);
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
		int limit = MakeNeoDB.maxFlankLength * 3;
		
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
				if(length >= limit) {
					paths.pollFirst();
					isRemoved = true;
				} else if(newExon.type == InputConvertorConstants.INS) {
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
						length -= path.refNucleotide.length();
					} else {
						length -= path.altNucleotide.length();
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
					length += newExon.refNucleotide.length();
				} else {
					length += newExon.altNucleotide.length();
				}
			}
			
			// reaching the start of transcript or exceeding length limitation
			if(length >= limit || newExon.start == -1) {
				ArrayList<Exon> mExons = new ArrayList<Exon>();
				Iterator<Exon> sExons = (Iterator<Exon>) paths.iterator();
				boolean isLeftMost = true;
				while(sExons.hasNext()) {
					Exon e = sExons.next();
					if(isLeftMost && newExon.start != -1) {
						e = e.copyExon(e.start + (length - limit), e.end);
						
						isLeftMost = false;
					}
					mExons.add(e);
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
					entry.sequence = leftSeq+mutExon.altNucleotide+rightSeq;
					entry.mutationMark = InputConvertorConstants.MUTATION_HEADER_ID;
					mutEntries.add(entry);
				}
			}
 		}
		
		return mutEntries;
	}
	
	
	
	public static String getSequenceFromExonList (List<Exon> exons) {
		StringBuilder sequence = new StringBuilder();
		
		for(Exon exon : exons) {
			if(exon.start == -1 || exon.start == Integer.MAX_VALUE) continue;
			
			if(exon.type == InputConvertorConstants.WILD) {
				sequence.append(exon.refNucleotide);
			} else {
				sequence.append(exon.altNucleotide);
			}
		}
		
		return sequence.toString();
	}
	
	public static ArrayList<FastaEntry> enumerateFastaEntry (GenomeLoader refGenome, Transcript t, Collection<Exon> exons, int times) {
		Exon[] exonGraph = Exon.toExonGraph(exons);
		ArrayList<FastaEntry> entries = new ArrayList<FastaEntry>();
		// if fail to generate sequence
		if(!refGenome.setSequence(t.chr, exonGraph)) {
			return entries;
		}
		
		ArrayList<FastaEntry> ntEntries = enumerateCombinatorialMutations(exonGraph);
		
		
		FastaEntry refEntry = new FastaEntry();
		
		StringBuilder seq = new StringBuilder();
		Exon exon = exonGraph[0];
		
		while(exon.start != Integer.MAX_VALUE) {
			if(exon.refNucleotide != null) {
				seq.append(exon.refNucleotide);
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
			
			for(int i=1; i<times; i++) {
				fastaEntry.sequence += fastaEntry.sequence;
				fastaEntry.description += fastaEntry.description;
			}
			
			for(int frame=0; frame<3; frame++) {
				FastaEntry aaEntry = new FastaEntry();
				String peptide = null;
				String tNucleotide = null;
				if(t.strand.equalsIgnoreCase("+")) {
					String[] seqs = Translator.translation(fastaEntry.sequence, frame);
					tNucleotide = seqs[0];
					peptide = seqs[1];
				} else if(t.strand.equalsIgnoreCase("-")) {
					String[] seqs = Translator.reverseComplementTranslation(fastaEntry.sequence, frame);
					tNucleotide = seqs[0];
					peptide = seqs[1];
				} else {
					// skip if there is no strand information
					continue;
				}
				aaEntry.description = fastaEntry.description;
				aaEntry.frame = frame+"";
				aaEntry.sequence = peptide;
				aaEntry.nucleotide = tNucleotide;
				aaEntry.transcript = t;
				aaEntry.mutationMark = fastaEntry.mutationMark;
				
				entries.add(aaEntry);
			}
			
		}
		
		return entries;
	}
	
	public static ArrayList<FastaEntry> enumerateFastaEntry (GenomeLoader refGenome, Transcript t, Collection<Exon> exons) {
		return enumerateFastaEntry(refGenome, t, exons, 1);
	}
	
	public static ArrayList<FastaEntry> enumerateFastaEntryCDS (GenomeLoader refGenome, Transcript t, Collection<Exon> exons) {
		Exon[] exonGraph = Exon.toExonGraph(exons);
		ArrayList<FastaEntry> entries = new ArrayList<FastaEntry>();
		// if fail to generate sequence
		if(!refGenome.setSequence(t.chr, exonGraph)) {
			return entries;
		}
		
		FastaEntry refEntry = new FastaEntry();
		
		StringBuilder seq = new StringBuilder();
		Exon exon = exonGraph[0];
		
		while(exon.start != Integer.MAX_VALUE) {
			if(exon.refNucleotide != null) {
				seq.append(exon.refNucleotide);
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
			String tNucleotide = null;
			if(t.strand.equalsIgnoreCase("+")) {
				String[] seqs = Translator.translation(refEntry.sequence, frame);
				tNucleotide = seqs[0];
				peptide = seqs[1];
			} else {
				String[] seqs = Translator.reverseComplementTranslation(refEntry.sequence, frame);
				tNucleotide = seqs[0];
				peptide = seqs[1];
			}
			
			String[] xPeptides = peptide.split("X");
			for(int i=0; i<xPeptides.length; i++) {
				if(longestLength < xPeptides[i].length()) {
					longestLength = xPeptides[i].length();
					longestFrameIdx = frame;
				}
			}
		}
		
		// write at most single mutated sequences
		getSingleMutation(exonGraph[0], 
				longestFrameIdx, t, 0,
				new LinkedList<Exon>(), entries);
		
		
		return entries;
	}
	
	private static void getSingleMutation (Exon exon, 
			int frame, Transcript t, int nCount,
			LinkedList<Exon> exonList, 
			ArrayList<FastaEntry> aaEntries) {
		
		if(exon.start == Integer.MAX_VALUE) {
			FastaEntry aaEntry = new FastaEntry();
			String nucleotide = getSequenceFromExonList(exonList);
			String peptide = null;
			String tNucleotide = null;
			if(t.strand.equalsIgnoreCase("+")) {
				String[] seqs = Translator.translation(nucleotide, frame);
				tNucleotide = seqs[0];
				peptide = seqs[1];
			} else {
				String[] seqs = Translator.reverseComplementTranslation(nucleotide, frame);
				tNucleotide = seqs[0];
				peptide = seqs[1];
			}
			aaEntry.description = getMutationDescription(exonList);
			aaEntry.frame = frame+"";
			aaEntry.sequence = peptide;
			aaEntry.nucleotide = tNucleotide;
			aaEntry.transcript = t;
			aaEntry.mutationMark = nCount > 0 ? 
					InputConvertorConstants.MUTATION_HEADER_ID : null;
			aaEntries.add(aaEntry);
			
			
		} else {
			for(Exon nExon : exon.nextExons) {
				if(nCount == 1) {
					// only wild type is acceptable.
					if(nExon.type == InputConvertorConstants.WILD) {
						exonList.add(nExon);
						getSingleMutation(nExon, frame, t, nCount, exonList, aaEntries);
						exonList.removeLast();
					}
				} else {
					exonList.add(nExon);
					if(nExon.type == InputConvertorConstants.WILD) {
						getSingleMutation(nExon, frame, t, nCount, exonList, aaEntries);
					} else {
						getSingleMutation(nExon, frame, t, nCount+1, exonList, aaEntries);
					}
					exonList.removeLast();
				}
			}
		}
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
	/*
	public static ArrayList<FastaEntry> removeDuplications (ArrayList<FastaEntry> entries) {
 		// it is possible to appear duplicated peptides because of several reference transcripts. 
 		Hashtable<String, FastaEntry> removeDuplication = new Hashtable<String, FastaEntry>();
 		ArrayList<FastaEntry> uniqueFastaEntries = new ArrayList<FastaEntry>();
 		int idx = 1;
         for(FastaEntry entry : entries) {
        	 FastaEntry firstEntry = removeDuplication.get(entry.getKey());
         	if(firstEntry != null) {
         		for(int i=0; i<entry.transcriptIds.size(); i++) {
         			String tID1 = entry.transcriptIds.get(i);
         			String tStrand1 = entry.strands.get(i);
         			
         			boolean isIncluded = false;
         			for(String tID2 : firstEntry.transcriptIds) {
         				if(tID1.equalsIgnoreCase(tID2)) {
         					isIncluded = true;
         				}
         			}
         			if(!isIncluded) {
         				firstEntry.transcriptIds.add(tID1);
         				firstEntry.strands.add(tStrand1);
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
	*/
}
