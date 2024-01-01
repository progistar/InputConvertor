package zhanglab.inputconvertor.data;

public class MS2PIPRecord implements Comparable<MS2PIPRecord>{

	public int idx;
	public String key;
	public String modifications;
	public String wildPeptide;
	public String charge;
	
	// for debug purpose.
	public String fullRecord;
	
	/**
	 * Only ic_peptide is supported.
	 * 
	 * @param icPeptide
	 * @return
	 */
	public static String getModifications (Peptide icPeptide) {
		/*********************************************************
		 * #RULE #HARD_CODING 
		 * IC_M_OXI to "m"
		 * 
		 * 
		 * 
		 *********************************************************/
		
		
		String modPeptide = icPeptide.modPeptide.replace(ModificationTable.IC_M_OXI, "m");
		String modifications = "";
		// I/L change
		// assume that the InferredPeptide has the same length to record.modifiedPeptide.
		for(int i=0; i<modPeptide.length(); i++) {
			char aa = modPeptide.charAt(i);
			
			// Oxidation
			if(aa == 'm') {
				if( modifications.length()!=0 ) {
					modifications += "|";
				}
				modifications += (i+1)+"|Oxidation";
			}
		}
		
		if(modifications.length() == 0) {
			modifications = "-";
		}
		
		// recover modification
		return modifications;
	}
	
	public static String getMS2PIPKey (String stripPeptide, String modifications, String charge) {
		return stripPeptide+"|"+modifications+"|"+charge;
	}

	@Override
	public int compareTo(MS2PIPRecord o) {
		
		if(this.idx < o.idx) {
			return -1;
		} else if(this.idx > o.idx) {
			return 1;
		}
		return 0;
	}
}
