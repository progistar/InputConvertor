package zhanglab.inputconvertor.data;

import java.util.Hashtable;

public class HarmonyData {

	// pipeline
	public String pipeline;
	public String label; // target and decoy
	public String spectrum;
	public String score;
	public String charge;
	
	
	// peptide information
	public String sequence;
	public String modification;
	public String aaVariant;
	
	// PSM key
	public String psmKey;
	
	// Genomic key
	public String genomicKey;
	
	// genomic annotation
	public String matchedLocation;
	public String matchedStrand;
	public String matchedMutation;
	public String matchedPeptide;
	public String matchedNucleotide;
	public String matchedReferenceNucleotide;
	public String matchedReadCount;
	public String matchedRPHM;
	public String proportion;
	public String geneId;
	public String geneName;
	public String geneStrand;
	public String geneType;
	public String classCode;
	public String uniqueClassCode;
	public String warningTag;
	
	// proteomic annotation
	public String proteinId;
	public String proteinGeneName;
	
	// penalty of class code
	public int penalty;
	public int proteinEvidence; // PE=1 is supposed to be "Reference sequence", otherwise non-reference.
	
	// fdr
	public String fdrClassId; //
	public String confidence; 
	
	// proteogenomic annotation
	public String majorGeneName;
	public String majorCategory;
	public String majorType;
	
	/**
	 * Major type and category:
	 * ** Type
	 * *** Reference
	 * IF or PE=1
	 * 
	 * *** Non-reference
	 * Else
	 * 
	 * 
	 * ** Category
	 * 
	 */
	
	public HarmonyData (String pipeline, String label, String spectrum, String score, String charge,
						String sequence, String modification,
						String matchedLocation, String matchedMutation, String matchedStrand,
						String matchedPeptide, String matchedNucleotide, String matchedReferenceNucleotide,
						String matchedReadCount, String matchedRPHM, String proportion,
						String geneId, String geneName, String geneStrand, String geneType,
						String classCode, String uniqueClassCode, String warningTag,
						String fdrClassId) {
		this.pipeline = pipeline;
		// PSM information
		this.label = label;
		this.spectrum = spectrum;
		this.sequence = sequence;
		this.modification = modification;
		this.score = score;
		this.charge = charge;
		// RNA-based sequence information
		this.matchedLocation = matchedLocation;
		this.matchedStrand = matchedStrand;
		this.matchedMutation = matchedMutation;
		this.matchedPeptide = matchedPeptide;
		this.matchedNucleotide = matchedNucleotide;
		this.matchedReferenceNucleotide = matchedReferenceNucleotide;
		// RNA quantification
		this.matchedReadCount = matchedReadCount;
		this.matchedRPHM = matchedRPHM;
		this.proportion = proportion;
		// Gene information
		this.geneId = geneId;
		this.geneName = geneName;
		this.geneStrand = geneStrand;
		this.geneType = geneType;
		// Event information
		this.classCode = classCode;
		this.uniqueClassCode = uniqueClassCode;
		this.warningTag = warningTag;
		// FDR class information
		this.fdrClassId = fdrClassId;
		
		
		// make a PSM key
		this.psmKey = this.label+"@"+this.spectrum+"@"+this.sequence+"@"+this.modification+"@"+this.charge;
		this.genomicKey = this.matchedLocation+"@"+this.matchedStrand+"@"+this.matchedMutation+"@"+this.matchedPeptide+"@"+this.matchedNucleotide;
		this.penalty = getPenalty(uniqueClassCode);
	}
	
	public String toString() {
		StringBuilder str = new StringBuilder();
		
		str.append(this.pipeline).append("\t")
		.append(this.spectrum).append("\t")
		.append(this.matchedPeptide).append("\t")
		.append(this.modification).append("\t")
		.append(this.charge).append("\t")
		.append(this.matchedPeptide.length()).append("\t")
		.append(this.score).append("\t")
		.append(this.confidence).append("\t")
		
		.append(this.matchedLocation).append("\t")
		.append(this.matchedStrand).append("\t")
		.append(this.matchedMutation).append("\t")
		.append(this.matchedNucleotide).append("\t")
		.append(this.matchedReferenceNucleotide).append("\t")
		.append(this.matchedReadCount).append("\t")
		.append(this.matchedRPHM).append("\t")
		.append(this.proportion).append("\t")
		
		.append(this.geneId).append("\t")
		.append(this.geneName).append("\t")
		.append(this.geneStrand).append("\t")
		.append(this.geneType).append("\t")
		.append(this.classCode).append("\t")
		.append(this.uniqueClassCode).append("\t")
		.append(this.warningTag).append("\t")
		
		.append(this.proteinId).append("\t")
		.append(this.proteinGeneName).append("\t")
		
		.append(this.majorGeneName).append("\t")
		.append(this.majorType).append("\t")
		.append(this.majorCategory);
		
		return str.toString();
	}
	

	private int getPenalty (String uniqueClassCode) {
		int penalty = 0;
		
		if(uniqueClassCode.contains("ES")) {
			penalty += 15;
		}
		
		if(uniqueClassCode.contains("ASS")) {
			penalty += 15;
		}
		
		if(uniqueClassCode.contains("UTR") || uniqueClassCode.contains("OOF")) {
			penalty += 20;
		}
		
		if(uniqueClassCode.contains("ncRNA")) {
			penalty += 30;
		}
		
		if(uniqueClassCode.contains("IR")) {
			penalty += 60;
		}
		
		if(uniqueClassCode.contains("asRNA")) {
			penalty += 120;
		}
		
		if(uniqueClassCode.contains("IGR")) {
			penalty += 240;
		}
		
		if(uniqueClassCode.contains("Unknown")) {
			penalty += 480;
		}
		
		return penalty;
	}
	
	public void setMajorAnnotation () {
		
		String proteomicType = "Reference";
		String genomicType = "Reference";
		
		if(this.proteinEvidence > 1) {
			proteomicType = "Non-reference";
		}
		if(!this.uniqueClassCode.equalsIgnoreCase("IF")) {
			genomicType = "Non-reference";
		}
		// point mutation (non-synonymous?)
		if(!this.matchedMutation.equals(".") && this.proteinEvidence > 1) {
			genomicType = "Non-reference";
		}
		
		if(genomicType.equalsIgnoreCase("Reference")) {
			this.majorGeneName = this.geneName;
			this.majorCategory = this.uniqueClassCode;
			this.majorType = genomicType;
		} 
		
		else if(proteomicType.equalsIgnoreCase("Reference")) {
			this.majorGeneName = this.proteinGeneName;
			this.majorCategory = "IF";
			this.majorType = proteomicType;
		}
		else {
			this.majorGeneName = this.geneName;
			this.majorCategory = this.uniqueClassCode;
			this.majorType = genomicType;
		}
		
		// remove duplicated major gene name
		String[] majorGeneNames = this.majorGeneName.split("\\|");
		this.majorGeneName = "";
		Hashtable<String, Boolean> check = new Hashtable<String, Boolean>();
		for(int i=0; i<majorGeneNames.length; i++) {
			if(check.get(majorGeneNames[i]) == null) {
				
				if(this.majorGeneName.length() > 0) {
					this.majorGeneName += "|";
				}
				this.majorGeneName += majorGeneNames[i];
				check.put(majorGeneNames[i], true);
			}
		}
	}
}
