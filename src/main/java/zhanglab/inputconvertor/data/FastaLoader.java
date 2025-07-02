package zhanglab.inputconvertor.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Hashtable;

import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;
import org.apache.commons.lang3.StringUtils;

import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.run.RunRefineNeoDB;

public class FastaLoader {

	public ArrayList<FastaEntry> entries = new ArrayList<FastaEntry>();
	
	public FastaLoader (File file) throws IOException {
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		FastaEntry entry = null;
		
		StringBuilder sequence = new StringBuilder();
		while((line = BR.readLine()) != null) {
			if(line.startsWith(">")) {
				
				if(entry != null) {
					entry.sequence = sequence.toString();
					sequence.setLength(0);
				}
				
				entry = new FastaEntry();
				entry.tool = InputConvertorConstants.REF_HEADER_ID;
				entry.idx = entries.size()+1;
				entry.originHeader = line.substring(1);
				entries.add(entry);
			} else {
				sequence.append(line);
			}
		}
		if(entry != null) {
			entry.sequence = sequence.toString();
		}
		
		BR.close();
	}
	
	public void removeSequenceOverlapped (ArrayList<FastaEntry> thatEntries, int partition) {
		ArrayList<String> sequences = new ArrayList<String>();
		Hashtable<String, String> removeList = new Hashtable<String, String>();
		// partition
		int junkSize = entries.size() / partition;
		for(int i=0; i<partition; i++) {
			int start = i * junkSize;
			int end = (i+1) * junkSize;
			if(i == partition -1) {
				end = entries.size();
			}
			
			
			for(int j=start; j<end; j++) {
				String[] peptides = StringUtils.split(entries.get(j).sequence.replace("I", "L"), AminoAcid.STOP_CODON_CHAR);
				for(String peptide : peptides) {
					if(removeList.get(peptide) == null) {
						sequences.add(peptide);
					}
				}
			}
			System.out.println((i+1)+" iteration ["+start+","+end+"): reading partial sequence "+sequences.size());
			Trie trie = Trie.builder().addKeywords(sequences).build();
			
			System.out.println("Done to build keyword-trie");
			
			for(FastaEntry refEntry : thatEntries) {
				Collection<Emit> matches = trie.parseText(refEntry.sequence.replace("I", "L"));
				for(Emit match : matches) {
					String overlap = match.getKeyword();
					removeList.put(overlap, "");
				}
			}
			
			sequences.clear();
		}
		
		
		ArrayList<FastaEntry> passEntries = new ArrayList<FastaEntry>();
		Hashtable<String, Boolean> compressEntires = new Hashtable<String, Boolean>();
		String[] headerMark = {"neodb"};
		for(FastaEntry entry : entries) {
			headerMark[0] = entry.originHeader.split("\\|")[0];
			
			String[] peptides = StringUtils.split(entry.sequence, AminoAcid.STOP_CODON_CHAR);
			boolean isPass = false;
			for(String peptide : peptides) {
				if(peptide.length() < InputConvertorConstants.MIN_PEPT_LEN) {
					continue;
				}
				
				if(removeList.get(peptide.replace("I", "L")) == null) {
					isPass = true;
					if(RunRefineNeoDB.isCompress) {
						compressEntires.put(peptide, true);
					}
				}
			}
			
			if(!RunRefineNeoDB.isCompress && isPass) {
				passEntries.add(entry);
			}
		}
		
		// if compress mode is on
		if(RunRefineNeoDB.isCompress) {
			int[] count = new int[1];
			compressEntires.forEach((peptide, nil)->{
				FastaEntry entry = new FastaEntry();
				entry.originHeader = headerMark[0]+"|"+headerMark[0]+"#"+(++count[0]);
				entry.sequence = peptide;
				passEntries.add(entry);
			});
		}
		
		if(RunRefineNeoDB.isCompress) {
			System.out.println("Compress: "+passEntries.size()); 
		} else {
			System.out.println("A total of "+(this.entries.size() - passEntries.size())+" were removed from the original entries");
			System.out.println("Remain: "+passEntries.size());
		}
		this.entries = passEntries;
	}
}
