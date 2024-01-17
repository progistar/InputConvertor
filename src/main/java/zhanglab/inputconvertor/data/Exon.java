package zhanglab.inputconvertor.data;

import java.util.ArrayList;

public class Exon implements Comparable<Exon> {

	public int start;
	public int end;
	public String nucleotide;
	public ArrayList<Mutation> snps = null;
	public ArrayList<Mutation> inss = null;
	public ArrayList<Mutation> dels = null;
	
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
	
	
}
