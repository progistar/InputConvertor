package zhanglab.inputconvertor.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Hashtable;

import zhanglab.inputconvertor.env.InputConvertorConstants;

public class GenomeLoader {

	private Hashtable<String, StringBuilder> genomeMap = new Hashtable<String, StringBuilder>();
	private VEPLoader vep = null;
	
	public GenomeLoader (File file) {
		// load genome
		try {
			
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String line = null;
			
			StringBuilder sequence = null;
			while((line = BR.readLine()) != null) {
				if(line.startsWith(">")) {
					String chr = line.split("\\s")[0].substring(1);
					System.out.println("Read "+chr);
					sequence = new StringBuilder();
					genomeMap.put(chr, sequence);
				} else {
					sequence.append(line.toUpperCase());
				}
			}
			
			BR.close();
			
		}catch(IOException ioe) {
			
		}
	}
	
	public void enrollVEPLaoder (VEPLoader vep) {
		System.out.println("## GenomeLoader ##");
		System.out.println("Enroll VEP");
		this.vep = vep;
	}
	
	/**
	 * zero-based, [start, end)
	 * 
	 * @param chr
	 * @param start
	 * @param end
	 * @return
	 */
	public String getSequence(String chr, int start, int end) {
		StringBuilder sequence = genomeMap.get(chr);
		if(sequence == null) {
			return null;
		}
		
		return sequence.substring(start, end);
	}
	
	public void setSequence(String chr, Collection<Exon> exons) {
		StringBuilder sequence = genomeMap.get(chr);
		for(Exon exon :exons) {
			int start = exon.start - 1;
			int end = exon.end;
			exon.nucleotide = sequence.subSequence(start,  end).toString();
			
			// enroll mutations
			if(vep != null) {
				ArrayList<Mutation> mutations = vep.getSNPByRange(chr, exon.start, exon.end+1);
				if(mutations.size() != 0) {
					exon.snps = mutations;
				}
				
				mutations = vep.getINSByRange(chr, exon.start, exon.end+1);
				if(mutations.size() != 0) {
					exon.inss = mutations;
				}
				
				mutations = vep.getDELByRange(chr, exon.start, exon.end+1);
				if(mutations.size() != 0) {
					exon.dels = mutations;
				}
			}
		}
	}
}
