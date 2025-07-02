package zhanglab.inputconvertor.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.input.HarmonyResult;

public class HarmonyFragPipe extends HarmonyResult {

	public HarmonyFragPipe(File file) throws IOException {
		super(file);
		this.pipeline = InputConvertorConstants.FRAGPIPE;
		renameHeader("[Score]", "Probability");
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		String[] fields = BR.readLine().split("\t");
		parseIndex(fields); // find pipeline.
		
		// parsing
		while((line = BR.readLine()) != null) {
			fields = line.split("\t");
			
			// PSM information
			String label = fields[labelIdx];
			// HARD-CODING for decoy prefix.
			if(label.startsWith("rev_")) {
				label = "-1";
			} else {
				label = "1";
			}
			
			String spectrum = fields[spectrumIdx];
			String sequence = fields[sequenceIdx];
			String modification = fields[modificationIdx];
			String score = fields[scoreIdx];
			String charge = fields[chargeIdx];
			
			// RNA-based sequence information
			String matchedLocation = fields[matchedLocationIdx];
			String matchedStrand = fields[matchedStrandIdx];
			String matchedMutation = fields[matchedMutationIdx];
			String matchedPeptide = fields[matchedPeptideIdx];
			String matchedNucleotide = fields[matchedNucleotideIdx];
			String matchedReferenceNucleotide = fields[matchedReferenceNucleotideIdx];
			// RNA quantification
			String matchedReadCount = fields[matchedReadCountIdx];
			String matchedRPHM = fields[matchedRPHMIdx];
			String proportion = fields[proportionIdx];
			// Gene information
			String geneId = fields[geneIdIdx];
			String geneName = fields[geneNameIdx];
			String geneStrand = fields[geneStrandIdx];
			String geneType = fields[geneTypeIdx];
			// Event information
			String classCode = fields[classCodeIdx];
			String uniqueClassCode = fields[uniqueClassCodeIdx];
			String warningTag = fields[warningTagIdx];		
			// FDR class id
			String fdrClassId = fields[fdrClassIdIdx];
			
			HarmonyData hData = new HarmonyData(this.pipeline, label, spectrum, score, charge, 
					sequence, modification, matchedLocation, matchedMutation, 
					matchedStrand, matchedPeptide, matchedNucleotide, matchedReferenceNucleotide, 
					matchedReadCount, matchedRPHM, proportion, 
					geneId, geneName, geneStrand, geneType, 
					classCode, uniqueClassCode, warningTag, fdrClassId);
			this.data.add(hData);
		}
		
		BR.close();
	}

	
	@Override
	public void parseIndex(String[] fields) {
		// TODO Auto-generated method stub
		super.parseIndex(fields);
		
		for(int i=0; i<fields.length; i++) {
			if(fields[i].equalsIgnoreCase("Spectrum") ) {
				this.spectrumIdx = i;
			} else if(fields[i].equalsIgnoreCase("Peptide") ) {
				this.sequenceIdx = i;
			} else if(fields[i].equalsIgnoreCase("Assigned Modifications") ) {
				this.modificationIdx = i;
			} else if(fields[i].equalsIgnoreCase("Probability") ) {
				this.scoreIdx = i;
			} else if(fields[i].equalsIgnoreCase("Charge") ) {
				this.chargeIdx = i;
			} else if(fields[i].equalsIgnoreCase("Protein") ) {
				this.labelIdx = i; // infer target decoy by "rev_" tag
			} else if(fields[i].equalsIgnoreCase("Class") ) {
				this.fdrClassIdIdx = i;
			}
		}
	}
}
