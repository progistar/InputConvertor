package zhanglab.inputconvertor.data;

import java.util.Hashtable;

public class Peptide {

	public String modPeptide;
	public String stripPeptide;
	public boolean isPass = true;
	
	/**
	 * Unified mod format:
	 * 
	 * AA + Mod-mass
	 * ex> C+57.021
	 * 
	 * 
	 * 
	 * @param peptide
	 * @param peptideFrom
	 */
	public Peptide (String peptide, Hashtable<String, String> ptmTable) {
		this.modPeptide = peptide;
		
		ptmTable.forEach((ptm, unimod)->{
			this.modPeptide = this.modPeptide.replace(ptm, unimod);
		});
		this.modPeptide = this.modPeptide.replace("-", "");
		this.stripPeptide = this.modPeptide.replaceAll("\\[(.*?)\\]","");
		
		// check if there are undesirable characters.
		for(int i=0; i<this.stripPeptide.length(); i++) {
			if(this.stripPeptide.charAt(i) - 'A' < 0 || this.stripPeptide.charAt(i) - 'A' > 26) {
				isPass = false;
			}
		}
	}
	
}
