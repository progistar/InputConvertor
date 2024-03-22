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
			this.modPeptide = peptide.replace(ModificationTable.PNOVO3_M_OXI, ModificationTable.IC_M_OXI);
			this.modPeptide = this.modPeptide.replace(ModificationTable.PNOVO3_C_CARBAM, ModificationTable.IC_C_CARBAM);
		}
		
		else if(peptideFrom.equalsIgnoreCase(InputConvertorConstants.PEAKS)) {
			this.modPeptide = peptide.replace(ModificationTable.PEAKS_M_OXI, ModificationTable.IC_M_OXI);
			this.modPeptide = this.modPeptide.replace(ModificationTable.PEAKS_C_CARBAM, ModificationTable.IC_C_CARBAM);
		}
		
		else if(peptideFrom.equalsIgnoreCase(InputConvertorConstants.COMET)) {
			this.modPeptide = peptide.replace(ModificationTable.COMET_M_OXI, ModificationTable.IC_M_OXI);
			this.modPeptide = this.modPeptide.replace(ModificationTable.COMET_C_CARBAM, ModificationTable.IC_C_CARBAM);
		}
		
		else if(peptideFrom.equalsIgnoreCase(InputConvertorConstants.AUTORT)) {
			this.modPeptide = peptide.replace(ModificationTable.AUTORT_M_OXI, ModificationTable.IC_M_OXI);
			this.modPeptide = this.modPeptide.replace(ModificationTable.AUTORT_C_CARBAM, ModificationTable.IC_C_CARBAM);
		}
		
		else if(peptideFrom.equalsIgnoreCase(InputConvertorConstants.PROSIT)) {
			this.modPeptide = peptide.replace(ModificationTable.UNIMOD_M_OXI, ModificationTable.IC_M_OXI);
			this.modPeptide = this.modPeptide.replace(ModificationTable.UNIMOD_C_CARBAM, ModificationTable.IC_C_CARBAM);
		}
		
		getStripFromModPeptide();
		setModifications();
	}
	
	private void setModifications () {
		String tempModPeptide = getLowerCaseMod(this.modPeptide);
		for(int i=0; i<tempModPeptide.length(); i++) {
			char aa = tempModPeptide.charAt(i);
			putMod(aa, i);
		}
	}
	
	public String getLowerCaseMod (String peptide) {
		String tempModPeptide = peptide.replace(ModificationTable.IC_M_OXI, "m");
		tempModPeptide = tempModPeptide.replace(ModificationTable.IC_C_CARBAM, "c");
		return tempModPeptide;
	}
	
	public String getUpperCaseMod (String peptide) {
		String tempModPeptide = peptide.replace("m", ModificationTable.IC_M_OXI);
		tempModPeptide = tempModPeptide.replace("c", ModificationTable.IC_C_CARBAM);
		return tempModPeptide;
	}
	
	private void putMod (char aa, int idx) {
		double modMass = 0;
		if(aa == 'm') {
			modMass = ModificationTable.OXIDATION_MASS;
		} else if(aa == 'c') {
			modMass = ModificationTable.CARBAM_MASS;
		}
		
		if(modMass != 0) {
			ArrayList<Double> mods = modifications.get(idx+1);
			if(mods == null) {
				mods = new ArrayList<Double>();
				modifications.put(idx+1, mods);
			}
			mods.add(modMass);
		}
	}
	
	private void getStripFromModPeptide () {
		this.stripPeptide = this.modPeptide.replace(ModificationTable.IC_M_OXI, "M");
		this.stripPeptide = this.stripPeptide.replace(ModificationTable.IC_C_CARBAM, "C");
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
		StringBuilder tempModPeptide = new StringBuilder(getLowerCaseMod(this.modPeptide));
		
		
		// I/L change
		// assume that the InferredPeptide has the same length to record.modifiedPeptide.
		for(int i=0; i<inferredPeptide.length(); i++) {
			char aa = inferredPeptide.charAt(i);
			if(aa == 'I' || aa == 'L') {
				tempModPeptide.setCharAt(i, aa);
			}
		}
		
		// restore both modPeptide and stripPeptide.
		this.modPeptide = getUpperCaseMod(tempModPeptide.toString());
		getStripFromModPeptide();
	}
	
	public void toAutoRTModPeptide () {
		this.modPeptide = this.modPeptide.replace(ModificationTable.IC_M_OXI, ModificationTable.AUTORT_M_OXI);
		this.modPeptide = this.modPeptide.replace(ModificationTable.IC_C_CARBAM, ModificationTable.AUTORT_C_CARBAM);
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
