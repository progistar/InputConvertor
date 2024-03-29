package zhanglab.inputconvertor.data;

import java.util.ArrayList;
import java.util.Collection;

import zhanglab.inputconvertor.env.InputConvertorConstants;

public class Exon implements Comparable<Exon> {

	public String chr;
	public int start;
	public int end;
	public String refNucleotide;
	public String altNucleotide;
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
	
	/**
	 * copy all except for next/prev connections
	 * @return
	 */
	public Exon copyExon () {
		Exon exon = new Exon(chr, start, end);
		exon.refNucleotide = this.refNucleotide;
		exon.altNucleotide = this.altNucleotide;
		exon.type = this.type;
		return exon;
	}
	
	public String getMutationDescription () {
		StringBuilder desc = new StringBuilder();
		desc.append(this.chr+":"+this.start+"-"+this.end);
		if(type != InputConvertorConstants.WILD) {
			
			if(type == InputConvertorConstants.SNP) {
				desc.append("[SNP]").append(this.refNucleotide).append(">").append(this.altNucleotide);
			} else if(type == InputConvertorConstants.MNP) {
				desc.append("[MNP]").append(this.refNucleotide).append(">").append(this.altNucleotide);
			} else if(type == InputConvertorConstants.DEL) {
				desc.append("[DEL]").append(this.refNucleotide);
			} else if(type == InputConvertorConstants.INS) {
				desc.append("[INS]").append(this.altNucleotide);
			}
			
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
	public static Exon[] toExonGraph (Collection<Exon> exons) {
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
		
		Exon[] returnExons = {startExon, endExon};
		
		return returnExons;
	}
	
	/**
	 * Divide exon only for wildtype
	 * 
	 * @param exon
	 * @param pos
	 * @return
	 */
	public static Exon divdeExon (Exon exon, int pos) {
		assert exon.type == InputConvertorConstants.WILD;
		
		ArrayList<Exon> prevs = exon.prevExons;
		ArrayList<Exon> nexts = exon.nextExons;
		
		
		// single size exon
		// or out of range
		if( (exon.end == pos) ||
			(exon.start > pos || exon.end < pos)) {
			return exon;
		}
		
		Exon leftExon = new Exon(exon.chr, exon.start, pos);
		leftExon.refNucleotide = exon.refNucleotide.substring(0, pos - exon.start + 1);
		
		Exon rightExon = new Exon(exon.chr, pos+1, exon.end);
		rightExon.refNucleotide = exon.refNucleotide.substring(pos - exon.start + 1, exon.end - exon.start + 1);
		
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
	
	public static ArrayList<Exon> mergeAdjacentExons (ArrayList<Exon> exons) {
		// type count
		boolean isUniqueType = true;
		for(int i=1; i<exons.size(); i++) {
			if(exons.get(i-1).type != exons.get(i).type) {
				isUniqueType = false;
			}
		}
		
		// only single-consistent type is allowed.
		assert isUniqueType;
		assert exons.size() != 0;
		
		ArrayList<Exon> adjExons = new ArrayList<Exon>();
		
		Exon firstExon = exons.get(0);
		Exon mExon = firstExon.copyExon();
		adjExons.add(mExon);
		
		for(int i=1; i<exons.size(); i++) {
			Exon exon = exons.get(i);
			// adjacent?
			if(exon.start == mExon.end+1) {
				mExon.end = exon.end;
				mExon.altNucleotide += exon.altNucleotide;
				mExon.refNucleotide += exon.refNucleotide;
			} else {
				mExon = exon.copyExon();
				adjExons.add(mExon);
			}
		}
		
		// make a connection
		for(int i=1; i<adjExons.size(); i++) {
			adjExons.get(i-1).nextExons.add(adjExons.get(i));
			adjExons.get(i).prevExons.add(adjExons.get(i-1));
		}
		
		return adjExons;
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
	
	public static Exon getWildTypeExon (Collection<Exon> exons) {
		Exon wildExon = null;
		
		for(Exon exon : exons) {
			if(exon.type == InputConvertorConstants.WILD) {
				wildExon = exon;
				break;
			}
		}
		
		return wildExon;
	}
}
