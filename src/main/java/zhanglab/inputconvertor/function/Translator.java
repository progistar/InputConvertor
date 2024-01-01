package zhanglab.inputconvertor.function;

import zhanglab.inputconvertor.data.Codon;

public class Translator {
	
	public static String translation (String nucleotides, int frame) {
		StringBuilder peptides = new StringBuilder();
		int length = nucleotides.length();
		for(int position=frame; position<length-2; position+=3) {
			char aa = Codon.nuclToAmino(nucleotides.substring(position,position+3));
			peptides.append(aa);
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
		for(int position=frame; position<nucleotides.length()-2; position+=3) {
			char aa = Codon.nuclToAmino(nucleotides.substring(position,position+3));
			peptides.append(aa);
		}
		return peptides.toString();
	}
}
