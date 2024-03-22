package zhanglab.inputconvertor.data;

import java.util.ArrayList;

public class Exon implements Comparable<Exon> {

	public String chr;
	public int start;
	public int end;
	public String nucleotide;
	public ArrayList<Mutation> snps = new ArrayList<Mutation>();
	public ArrayList<Mutation> inss = new ArrayList<Mutation>();;
	public ArrayList<Mutation> dels = new ArrayList<Mutation>();
	
	public Exon (String chr, int start, int end) {
		this.chr = chr;
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
				int relPos = snp.pos - this.start;
				mutantSeq.setCharAt(relPos, snp.altSeq.charAt(0));
			}
		}
		
		return mutantSeq.toString();
	}
	
	public String getSequenceWithINS (String insID, boolean germSNP, boolean somaSNP) {
		String sequence = getSequenceWithSNPs(germSNP, somaSNP);
		Mutation ins = null;
		for(Mutation mutation : inss) {
			if(mutation.key.equalsIgnoreCase(insID)) {
				ins = mutation;
			}
		}
		
		if(ins == null) {
			return sequence;
		}
		
		StringBuilder mutantSeq = new StringBuilder();
		
		for(int i=start; i<=end; i++) {
			int relPos = i-start;
			mutantSeq.append(sequence.charAt(relPos));
			if(i == ins.pos) {
				mutantSeq.append(ins.altSeq);
			}
		}
		
		return mutantSeq.toString();
	}
	
	public String getSequenceWithDEL (String delID, boolean germSNP, boolean somaSNP) {
		String sequence = getSequenceWithSNPs(germSNP, somaSNP);
		Mutation del = null;
		for(Mutation mutation : dels) {
			if(mutation.key.equalsIgnoreCase(delID)) {
				del = mutation;
			}
		}
		
		if(del == null) {
			return sequence;
		}
		
		StringBuilder mutantSeq = new StringBuilder();
		
		for(int i=start; i<=end; i++) {
			int relPos = i-start;
			if(i == del.pos) {
				mutantSeq.append(del.altSeq);
			} else {
				mutantSeq.append(sequence.charAt(relPos));
			}
		}
		
		
		return mutantSeq.toString();
	}
	
	public String getMutationDescription (String indelID, boolean germSNP, boolean somaSNP) {
		StringBuilder desc = new StringBuilder(this.chr+":"+this.start+"-"+this.end);
		
		desc.append("(");
		int mutationCnt = 0;
		for(Mutation snp : snps) {
			if((snp.isSomatic && somaSNP) || (!snp.isSomatic && germSNP)) {
				if(mutationCnt != 0) {
					desc.append(",");
				}
				desc.append(snp.key);
				mutationCnt++;
			}
		}
		
		if(indelID != null) {
			for(Mutation mutation : dels) {
				if(mutation.key.equalsIgnoreCase(indelID)) {
					if(mutationCnt != 0) {
						desc.append(",");
					}
					desc.append(mutation.key);
					mutationCnt++;
					break;
				}
			}
			
			for(Mutation mutation : inss) {
				if(mutation.key.equalsIgnoreCase(indelID)) {
					if(mutationCnt != 0) {
						desc.append(",");
					}
					desc.append(mutation.key);
					mutationCnt++;
					break;
				}
			}
		}
		
		if(mutationCnt == 0) {
			desc.append("-");
		}
		
		desc.append(")");
		return desc.toString();
	}
	
}
