package zhanglab.inputconvertor.function;

import java.util.Collection;

import zhanglab.inputconvertor.data.Codon;
import zhanglab.inputconvertor.data.Exon;
import zhanglab.inputconvertor.env.InputConvertorConstants;

public class Translator {
	
	public static String translation (Collection<Exon> exons, int frame) {
		StringBuilder nucleotide = new StringBuilder();
		for(Exon exon : exons) {
			nucleotide.append(exon.nucleotide);
		}
		
		return translation(nucleotide.toString(), frame);
	}
	
	public static String reverseComplementTranslation (Collection<Exon> exons, int frame) {
		StringBuilder nucleotide = new StringBuilder();
		for(Exon exon : exons) {
			nucleotide.append(exon.nucleotide);
		}
		return reverseComplementTranslation(nucleotide.toString(), frame);
	}
	
	public static String translation (String nucleotides, int frame) {
		StringBuilder peptides = new StringBuilder();
		int length = nucleotides.length();
		StringBuilder codon = new StringBuilder();
		for(int pos=frame; pos<length; pos++) {
			
			char nt = nucleotides.charAt(pos);
			
			// skip deletion
			if(nt == InputConvertorConstants.DELETION_MARK.charAt(0)) {
				continue;
			}
			
			codon.append(nt);
			if(codon.length() == 3) {
				String codonStr = codon.toString();
				char aa = Codon.nuclToAmino(codonStr.toUpperCase());
				peptides.append(aa);
				codon.setLength(0);
			}
		}
		return peptides.toString();
	}
	
	public static String reverseComplementTranslation (String nucleotides, int frame) {
		StringBuilder peptides = new StringBuilder();
		StringBuilder reverseComplementNTs = new StringBuilder(nucleotides);
		int length = nucleotides.length();
		for(int i=0; i<length; i++) {
			switch(reverseComplementNTs.charAt(i)) {
				case 'A': reverseComplementNTs.setCharAt(i, 'T'); break;
				case 'C': reverseComplementNTs.setCharAt(i, 'G'); break;
				case 'T': reverseComplementNTs.setCharAt(i, 'A'); break;
				case 'G': reverseComplementNTs.setCharAt(i, 'C'); break;
				default : break;
			}
		}
		
		nucleotides = reverseComplementNTs.reverse().toString();
		
		StringBuilder codon = new StringBuilder();
		for(int pos=frame; pos<length; pos++) {
			
			char nt = nucleotides.charAt(pos);
			
			// skip deletion
			if(nt == InputConvertorConstants.DELETION_MARK.charAt(0)) {
				continue;
			}
			
			codon.append(nt);
			if(codon.length() == 3) {
				String codonStr = codon.toString();
				char aa = Codon.nuclToAmino(codonStr.toUpperCase());
				peptides.append(aa);
				codon.setLength(0);
			}
		}
		
		return peptides.toString();
	}
}
