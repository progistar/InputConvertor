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

import zhanglab.inputconvertor.env.InputConvertorConstants;

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
				String[] peptides = entries.get(j).sequence.replace("I", "L").split("X");
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
		for(FastaEntry entry : entries) {
			String[] peptides = entry.sequence.replace("I", "L").split("X");
			boolean isPass = false;
			for(String peptide : peptides) {
				if(peptide.length() == 0) {
					continue;
				}
				
				if(removeList.get(peptide) == null) {
					isPass = true;
					break;
				}
			}
			if(isPass) {
				passEntries.add(entry);
			}
		}
		
		System.out.println("A total of "+(this.entries.size() - passEntries.size())+" were removed from the original entries");
		System.out.println("Remain: "+passEntries.size());
		this.entries = passEntries;
	}
}
