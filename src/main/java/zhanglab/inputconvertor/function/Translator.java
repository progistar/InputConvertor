package zhanglab.inputconvertor.function;

import zhanglab.inputconvertor.data.Codon;
import zhanglab.inputconvertor.env.InputConvertorConstants;

public class Translator {
	
	public static String[] translation (String nucleotides, int frame) {
		StringBuilder peptides = new StringBuilder();
		int length = nucleotides.length();
		
		StringBuilder translatedNucleotides = new StringBuilder();
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
				translatedNucleotides.append(codonStr);
				codon.setLength(0);
			}
		}
		
		String[] newString = {translatedNucleotides.toString(), peptides.toString()};
		
		return newString;
	}
	
	public static String getReverseComplement (String nucleotides) {
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
		
		return reverseComplementNTs.reverse().toString();
	}
	
	public static String[] reverseComplementTranslation (String nucleotides, int frame) {
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
		
		reverseComplementNTs.setLength(0);
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
				reverseComplementNTs.append(codonStr);
				codon.setLength(0);
			}
		}
		
		String[] newString = {reverseComplementNTs.toString(), peptides.toString()};
		
		return newString;
	}
}
