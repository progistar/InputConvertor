package zhanglab.inputconvertor.data;

import java.util.ArrayList;

public class Exon implements Comparable<Exon> {

	public int start;
	public int end;
	public String nucleotide;
	public ArrayList<Mutation> snps = new ArrayList<Mutation>();
	public ArrayList<Mutation> inss = new ArrayList<Mutation>();;
	public ArrayList<Mutation> dels = new ArrayList<Mutation>();
	
	public Exon (int start, int end) {
		this.start = start;
		this.end = end;
	}

	@Override
	public int compareTo(Exon o) {
		if(this.start < o.start) {
			return -1;
		}else if(this.start > o.start) {
			return 1;
		}
		return 0;
	}
	
	public boolean isMutant () {
		boolean isMutant = false;
		
		if(snps.size() > 0 || inss.size() > 0 || dels.size() > 0) {
			isMutant = true;
		}
		
		return isMutant;
	}
	
	public String getSequenceWithSNPs (boolean germSNP, boolean somaSNP) {
		if(snps.size() == 0) return this.nucleotide;
		
		StringBuilder mutantSeq = new StringBuilder(this.nucleotide);
		
		for(Mutation snp : snps) {
			if((snp.isSomatic && somaSNP) || (!snp.isSomatic && germSNP)) {
				mutantSeq.setCharAt(snp.pos, snp.altSeq.charAt(0));
			}
		}
		
		return mutantSeq.toString();
	}
	
	public String getSequenceWithINS (String key) {
		Mutation ins = null;
		for(Mutation mutation : inss) {
			if(mutation.key.equalsIgnoreCase(key)) {
				ins = mutation;
			}
		}
		
		if(ins == null) {
			return this.nucleotide;
		}
		
		StringBuilder mutantSeq = new StringBuilder();
		
		for(int i=start; i<=end; i++) {
			int relPos = i-start;
			mutantSeq.append(this.nucleotide.charAt(relPos));
			if(i == ins.pos) {
				mutantSeq.append(ins.altSeq);
			}
		}
		
		
		return mutantSeq.toString();
	}
	
	public String getSequenceWithDEL (String key) {
		Mutation ins = null;
		for(Mutation mutation : inss) {
			if(mutation.key.equalsIgnoreCase(key)) {
				ins = mutation;
			}
		}
		
		if(ins == null) {
			return this.nucleotide;
		}
		
		StringBuilder mutantSeq = new StringBuilder();
		
		for(int i=start; i<=end; i++) {
			int relPos = i-start;
			mutantSeq.append(this.nucleotide.charAt(relPos));
			if(i == ins.pos) {
				mutantSeq.append(ins.altSeq);
			}
		}
		
		
		return mutantSeq.toString();
	}
	
}
