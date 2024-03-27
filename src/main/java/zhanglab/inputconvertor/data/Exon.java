package zhanglab.inputconvertor.data;

import java.util.ArrayList;
import java.util.Collection;

import zhanglab.inputconvertor.env.InputConvertorConstants;

public class Exon implements Comparable<Exon> {

	public String chr;
	public int start;
	public int end;
	public String nucleotide;
	public Mutation mutation;
	public byte type = InputConvertorConstants.WILD;
	public ArrayList<Exon> nextExons = new ArrayList<Exon>();
	public ArrayList<Exon> prevExons = new ArrayList<Exon>();
	
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
	
	public Exon copyExon () {
		Exon exon = new Exon(chr, start, end);
		exon.nucleotide = this.nucleotide;
		return exon;
	}
	
	public String getMutationDescription () {
		StringBuilder desc = new StringBuilder(this.chr+":"+this.start+"-"+this.end);
		if(this.mutation != null) {
			desc.append("(").append(this.mutation.key).append(")");
		}
		return desc.toString();
	}
	
	public boolean hasPosition (int pos) {
		if(this.start <= pos && this.end >= pos) {
			return true;
		}
		
		return false;
	}
	
	/**
	 * 
	 * Return new exon list
	 * 
	 * @param exons
	 * @return
	 */
	public static ArrayList<Exon> addStartEndEonxs (Collection<Exon> exons) {
		ArrayList<Exon> newExons = new ArrayList<Exon>();
		
		Exon startExon = new Exon("", -1, -1);
		Exon endExon = new Exon("", Integer.MAX_VALUE, Integer.MAX_VALUE);
		newExons.add(startExon);
		for(Exon exon : exons) {
			newExons.add(exon.copyExon());
			
		}
		newExons.add(endExon);
		

		for(int i=1; i<newExons.size(); i++) {
			newExons.get(i-1).nextExons.add(newExons.get(i));
			newExons.get(i).prevExons.add(newExons.get(i-1));
		}
		
		return newExons;
	}
	
	public static Exon divdeExon (Exon exon, int pos) {
		ArrayList<Exon> prevs = exon.prevExons;
		ArrayList<Exon> nexts = exon.nextExons;
		
		
		// single size exon
		// or out of range
		if( (exon.end == pos) ||
			(exon.start > pos || exon.end < pos)) {
			return exon;
		}
		
		Exon leftExon = new Exon(exon.chr, exon.start, pos);
		leftExon.nucleotide = exon.nucleotide.substring(0, pos - exon.start + 1);
		
		Exon rightExon = new Exon(exon.chr, pos+1, exon.end);
		rightExon.nucleotide = exon.nucleotide.substring(pos - exon.start + 1, exon.end - exon.start + 1);
		
		leftExon.nextExons.add(rightExon);
		rightExon.prevExons.add(leftExon);
		
		// disconnect and reconnect
		for(Exon pExon : prevs) {
			// disconnect to old exon
			disconnectTargetExon(pExon, exon);
			
			// connect to new exon
			pExon.nextExons.add(leftExon);
			leftExon.prevExons.add(pExon);
		}
		for(Exon nExon : nexts) {
			disconnectTargetExon(nExon, exon);
			
			// connect to new exon
			nExon.prevExons.add(rightExon);
			rightExon.nextExons.add(nExon);
		}
		
		return leftExon;
	}
	
	public static void disconnectTargetExon (Exon exon, Exon tExon) {
		int pCnt = 0;		int nCnt = 0; 
		while(exon.prevExons.remove(tExon)) {pCnt++;}
		while(exon.nextExons.remove(tExon)) {nCnt++;}
		if(pCnt > 1 || nCnt > 1) {
			System.out.println("Duplicated connection is detected");
			System.out.println("prve dups: "+pCnt+"& next dups: "+nCnt);
		}
	}
	
	public static boolean isConnectedList(Collection<Exon> exons) {
		boolean isConnected = false;
		
		for(Exon exon : exons) {
			if(exon.start == -1) {
				isConnected = true;
			}
			break;
		}
		
		return isConnected;
	}
}
