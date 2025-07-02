package zhanglab.inputconvertor.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.input.HarmonyResult;

public class HarmonypXg extends HarmonyResult {

	public HarmonypXg(File file) throws IOException {
		super(file);
		this.pipeline = InputConvertorConstants.PXG;
		renameHeader("[Score]", "Percolator_score");
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		// FINDA FOLLOWING PATTERNS: AGBD[UNIMOD:35]AMDK+358AMKS(+382)AMDK[+382]AADDK+358.3AMKS(+382.1)AMDK[+382.2]AAD
		String ptmParserRegExr			=	"(\\[\\w+:\\d+\\]|\\([+-]?\\d+\\.*\\d+\\)|[+-]?\\d+\\.*\\d+|\\[[+-]?\\d+\\.*\\d+\\])";
		Pattern ptmPattern = Pattern.compile(ptmParserRegExr);
		
		String[] fields = BR.readLine().split("\t");
		parseIndex(fields); // find pipeline.
		
		// parsing
		while((line = BR.readLine()) != null) {
			fields = line.split("\t");
			
			// PSM information
			String label = fields[labelIdx];
			String spectrum = fields[spectrumIdx];
			String sequence = fields[sequenceIdx];
			String aaVariant = fields[aaVariantIdx];
			String modification = "";
			
			// parse inferred peptide
			// check PTM pattern
			Matcher matcher = ptmPattern.matcher(sequence);
			while(matcher.find()) {
				if(modification.length() > 0) {
					modification += "|";
				}
				modification += matcher.group();
			}
			// add AAVariant
			// this is because length 1 should be "-" or "." (empty field).
			if(aaVariant.length() > 1) {
				if(modification.length() > 0) {
					modification += "|";
				}
				modification += aaVariant;
			}
			
			// remove unimod and mass relating patterns
			sequence = sequence.replaceAll(ptmParserRegExr, "");
			
			
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
					sequence, modification, 
					matchedLocation, matchedMutation, matchedStrand, 
					matchedPeptide, matchedNucleotide, matchedReferenceNucleotide, 
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
			if(pipeline.equalsIgnoreCase(InputConvertorConstants.PXG)) {
				if(fields[i].equalsIgnoreCase("SpecID") ) {
					this.spectrumIdx = i;
				} else if(fields[i].equalsIgnoreCase("InferredPeptide") ) {
					this.sequenceIdx = i;
				} else if(fields[i].equalsIgnoreCase("AminoAcidVariant") ) {
					this.aaVariantIdx = i;
				} else if(fields[i].equalsIgnoreCase("final_score") ) {
					this.scoreIdx = i;
				} else if(fields[i].equalsIgnoreCase("ic_charge") ) {
					this.chargeIdx = i;
				} else if(fields[i].equalsIgnoreCase("Label") ) {
					this.labelIdx = i;
				} else if(fields[i].equalsIgnoreCase("IsReference") ) {
					this.fdrClassIdIdx = i;
				}
			}
		}
	}
}
