package zhanglab.inputconvertor.data;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;

import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.env.ProteomeConstants;

public class Peptide {

	public String modPeptide;
	public String stripPeptide;
	public Hashtable<Integer, ArrayList<Double>> modifications =  new Hashtable<Integer, ArrayList<Double>>();
	
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
	public Peptide (String peptide, String peptideFrom) {
		if(peptideFrom.equalsIgnoreCase(InputConvertorConstants.CASANOVO) ||
			peptideFrom.equalsIgnoreCase(InputConvertorConstants.IC_CONVERTOR)) {
			this.modPeptide = peptide;
		}
		
		else if(peptideFrom.equalsIgnoreCase(InputConvertorConstants.PNOVO3)) {
			String M = ModificationTable.PNOVO3_M_OXI;
			this.modPeptide = peptide.replace(M, ModificationTable.IC_M_OXI);
		}
		
		else if(peptideFrom.equalsIgnoreCase(InputConvertorConstants.PEAKS)) {
			String M = ModificationTable.PEAKS_M_OXI;
			this.modPeptide = peptide.replace(M, ModificationTable.IC_M_OXI);
		}
		
		else if(peptideFrom.equalsIgnoreCase(InputConvertorConstants.COMET)) {
			String M = ModificationTable.COMET_M_OXI;
			this.modPeptide = peptide.replace(M, ModificationTable.IC_M_OXI);
			
			// select peptide part only
			// K.ACM+15.995GG.A => ACM+15.995GG
			int firstIdx = this.modPeptide.indexOf(".");
			int lastIdx = this.modPeptide.lastIndexOf(".");
			this.modPeptide = this.modPeptide.substring(firstIdx+1, lastIdx);
		}
		
		else if(peptideFrom.equalsIgnoreCase(InputConvertorConstants.AUTORT)) {
			String M = ModificationTable.AUTORT_M_OXI;
			this.modPeptide = peptide.replace(M, ModificationTable.IC_M_OXI);
			
		}
		
		getStripFromModPeptide();
		setModifications();
	}
	
	private void setModifications () {
		String tempModPeptide = this.modPeptide.replace(ModificationTable.IC_M_OXI, "m");
		
		for(int i=0; i<tempModPeptide.length(); i++) {
			char aa = tempModPeptide.charAt(i);
			if(aa == 'm') {
				ArrayList<Double> mods = modifications.get(i+1);
				if(mods == null) {
					mods = new ArrayList<Double>();
					modifications.put(i+1, mods);
				}
				mods.add(ModificationTable.OXIDATION_MASS);
			}
		}
	}
	
	private void getStripFromModPeptide () {
		this.stripPeptide = this.modPeptide.replace(ModificationTable.IC_M_OXI, "M");
	}
	
	public double getMass (boolean withMods) {
		double mass = ProteomeConstants.H2O;
		int sequenceLength = this.stripPeptide.length();
		for(int i=0; i<sequenceLength; i++) {
			
			mass += AminoAcid.getAminoAcid(this.stripPeptide.charAt(i)).getResidualMass();
			
			ArrayList<Double> mods = this.modifications.get(i+1);
			if(withMods && mods != null) { 
				for(Double modMass : mods) {
					mass += modMass;
				}
			}
		}
		
		return mass;
	}

	public void icPeptideToILDetermination (String inferredPeptide) {
		StringBuilder tempModPeptide = new StringBuilder(this.modPeptide.replace(ModificationTable.IC_M_OXI, "m"));
		
		// I/L change
		// assume that the InferredPeptide has the same length to record.modifiedPeptide.
		for(int i=0; i<inferredPeptide.length(); i++) {
			char aa = inferredPeptide.charAt(i);
			if(aa == 'I' || aa == 'L') {
				tempModPeptide.setCharAt(i, aa);
			}
		}
		
		// restore both modPeptide and stripPeptide.
		this.modPeptide = tempModPeptide.toString().replace("m", ModificationTable.IC_M_OXI);
		getStripFromModPeptide();
	}
	
	public void toAutoRTModPeptide () {
		this.modPeptide = this.modPeptide.replace(ModificationTable.IC_M_OXI, ModificationTable.AUTORT_M_OXI);
	}
	
	public double[] getTheoreticalLadder (double ion, double charge, boolean withMods) {
		int sequenceLength = this.stripPeptide.length();
		double[] ladder = new double[sequenceLength];
		
		double nTermModi = 0;
		for(int i=0; i<sequenceLength; i++) {
			ladder[i] = AminoAcid.getAminoAcid(this.stripPeptide.charAt(i)).getResidualMass();
			ArrayList<Double> mods = this.modifications.get(i+1);
			if(withMods && mods != null) { 
				for(Double modMass : mods) {
					ladder[i] += modMass;
				}
			}
		}
		
		// x, y, z ions
		if(ion == ProteomeConstants.X_ION || ion == ProteomeConstants.Y_ION || ion == ProteomeConstants.Z_ION) {
			// accumulated ions
			for(int i=sequenceLength-2; i>=0; i--) ladder[i] += ladder[i+1];
		}
		// a, b, c ions
		else if (ion == ProteomeConstants.A_ION || ion == ProteomeConstants.B_ION || ion == ProteomeConstants.C_ION){
			// accumulated ions
			ladder[0] += nTermModi;
			for(int i=1; i<sequenceLength; i++) ladder[i] += ladder[i-1];
		} else {
			// unsupported ions
			System.err.println("Unsupported IONs: "+ion);
			return null;
		}
		
		// charge and ion
		for(int i=0; i<sequenceLength; i++) ladder[i] = (ladder[i] + (charge-1)*ProteomeConstants.Proton + ion) / charge;
		
		Arrays.sort(ladder);
		
		double[] refinedLadder = new double[ladder.length-1];
		for(int i=0; i<refinedLadder.length; i++) {
			refinedLadder[i] = ladder[i];
		}
		
		return refinedLadder;
	}
}
