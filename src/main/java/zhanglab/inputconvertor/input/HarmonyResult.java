package zhanglab.inputconvertor.input;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;

import zhanglab.inputconvertor.data.FastaEntry;
import zhanglab.inputconvertor.data.FastaLoader;
import zhanglab.inputconvertor.data.HarmonyData;
import zhanglab.inputconvertor.data.HarmonyFragPipe;
import zhanglab.inputconvertor.data.HarmonyProtein;

abstract public class HarmonyResult {
	public ArrayList<HarmonyData> data = new ArrayList<HarmonyData>();

	// predefined header
	public static String[] HEADER = {
			"Pipeline",
			
			"Spectrum",
			"Sequence", // matchedPeptide
			"Modification", 
			"Charge",
			"Length",
			"[Score]",
			"Confidence",
			
			"Genomic_location",
			"Strand",
			"Mutation",
			"Observed_nucleotide",
			"Reference_nucleotide",
			"Read_count",
			"RPHM",
			"Proportion",
			
			"Gene_id",
			"Gene_name",
			"Gene_strand",
			"Gene_type",
			"Class_code",
			"Unique_class_code",
			"Warning_tag",
			
			"Protein_id",
			"Protein_gene_name",
			
			"Major_gene_name",
			"Major_type",
			"Major_category"
	};
	
	public String pipeline = "Undefined";
	// PSM
	public int labelIdx = -1;
	public int spectrumIdx = -1; // pXg: SpecID
	public int sequenceIdx = -1;
	public int scoreIdx = -1; // FragPipe: probability
	public int chargeIdx = -1;
	public int aaVariantIdx = -1;
	public int modificationIdx = -1;
	
	// genomic annotation index
	public int matchedLocationIdx = -1;
	public int matchedMutationIdx = -1;
	public int matchedStrandIdx = -1;
	public int matchedPeptideIdx = -1;
	public int matchedNucleotideIdx = -1;
	public int matchedReferenceNucleotideIdx = -1;
	public int matchedReadCountIdx = -1;
	public int matchedRPHMIdx = -1;
	public int proportionIdx = -1;
	public int geneIdIdx = -1;
	public int geneNameIdx = -1;
	public int geneStrandIdx = -1;
	public int geneTypeIdx = -1;
	public int classCodeIdx = -1;
	public int uniqueClassCodeIdx = -1;
	public int warningTagIdx = -1;
	
	// fdr class id
	public int fdrClassIdIdx = -1;
	
	
	// example
	public static void main(String[] args) throws IOException {
		HarmonyResult r = null;
		
		File file = new File("/Volumes/Papers/2025_NeoFlow2/test_proteogenomic_harmonization/C3L-04911_I.annotate.tsv");
		r = new HarmonyFragPipe(file);
		//File file = new File("/Volumes/Papers/2025_NeoFlow2/test_proteogenomic_harmonization/DN_C3L-04911_I.annotate.tsv");
		//r = new HarmonypXg(file);
		
		File ref = new File("/Volumes/Papers/2024_TNBC/3.Database/non_compressed/SUM159.pos.neodb.fasta");
		FastaLoader fLoader = new FastaLoader(ref);
		
		r.runSequentialSteps(fLoader, "rev_");
		
		for(int i=0; i<HEADER.length; i++) {
			System.out.print(HEADER[i]+"\t");
		}
		System.out.println();
		for(int i=0; i<r.data.size(); i++) {
			System.out.println(r.data.get(i).toString());
		}
	}
	
	public HarmonyResult (File file) throws IOException {
		System.out.println("Parse result: "+file.getName());
	}
	
	/**
	 * Check pipeline. <br>
	 * pXg or FragPipe.
	 * 
	 * @param fields
	 * @return
	 */
	public void parseIndex (String[] fields) {
		
		// infer pipeline using predefined fields.
		for(int i=0; i<fields.length; i++) {
			// genomic index
			if(fields[i].equalsIgnoreCase("Matched_location") ) {
				this.matchedLocationIdx = i;
			} else if(fields[i].equalsIgnoreCase("Matched_mutations") ) {
				this.matchedMutationIdx = i;
			} else if(fields[i].equalsIgnoreCase("Matched_strand") ) {
				this.matchedStrandIdx = i;
			} else if(fields[i].equalsIgnoreCase("Matched_peptide") ) {
				this.matchedPeptideIdx = i;
			} else if(fields[i].equalsIgnoreCase("Matched_nucleotide") ) {
				this.matchedNucleotideIdx = i;
			} else if(fields[i].equalsIgnoreCase("Matched_reference_nucleotide") ) {
				this.matchedReferenceNucleotideIdx = i;
			} else if(fields[i].equalsIgnoreCase("Matched_read_count") ) {
				this.matchedReadCountIdx = i;
			} else if(fields[i].equalsIgnoreCase("Matched_RPHM") ) {
				this.matchedRPHMIdx = i;
			} else if(fields[i].equalsIgnoreCase("Proportion") ) {
				this.proportionIdx = i;
			} else if(fields[i].equalsIgnoreCase("Gene_id") ) {
				this.geneIdIdx = i;
			} else if(fields[i].equalsIgnoreCase("Gene_name") ) {
				this.geneNameIdx = i;
			} else if(fields[i].equalsIgnoreCase("Gene_strand") ) {
				this.geneStrandIdx = i;
			} else if(fields[i].equalsIgnoreCase("Gene_type") ) {
				this.geneTypeIdx = i;
			} else if(fields[i].equalsIgnoreCase("Class_code") ) {
				this.classCodeIdx = i;
			} else if(fields[i].equalsIgnoreCase("Unique_class_code") ) {
				this.uniqueClassCodeIdx = i;
			} else if(fields[i].equalsIgnoreCase("Warning_tag") ) {
				this.warningTagIdx = i;
			}
		}
	}
	/**
	 * 
	 * 
	 * 
	 * @param fasta
	 * @param decoyPrefix
	 */
	public void runSequentialSteps (FastaLoader fasta, String decoyPrefix) {
		// prioritize major "genomic" form.
		this.prioritizeMajorGenomicForm();
		// assign protein id and get major annotation.
		this.assignProteinId(fasta, decoyPrefix);
		// assign confidence level
		this.assignConfidence();
		// drop decoy labeled data
		this.dropDecoyLabel();
	}
	
	/**
	 * Sequential step 1:
	 * 
	 */
	private void prioritizeMajorGenomicForm () {
		// 1. Sort by proportion and penalty
		Collections.sort(this.data, new Comparator<HarmonyData>() {
			@Override
			public int compare(HarmonyData o1, HarmonyData o2) {
				double p1 = Double.parseDouble(o1.proportion);
				double p2 = Double.parseDouble(o2.proportion);
				if(p1 > p2) {
					return -1;
				} else if(p1 < p2) {
					return 1;
				}
				return 0;
			}
		});
		
		Collections.sort(this.data, new Comparator<HarmonyData>() {
			@Override
			public int compare(HarmonyData o1, HarmonyData o2) {
				if(o1.penalty < o2.penalty) {
					return -1;
				} else if(o1.penalty > o2.penalty) {
					return 1;
				}
				return 0;
			}
		});
		
		// Select top ranked one
		ArrayList<HarmonyData> rankedData = new ArrayList<HarmonyData>();
		
		Hashtable<String, Boolean> check = new Hashtable<String, Boolean>();
		for(int i=0; i<this.data.size(); i++) {
			HarmonyData hData = this.data.get(i);
			
			if(check.get(hData.psmKey) == null) {
				rankedData.add(hData);
				check.put(hData.psmKey, true);
			}
		}
		
		this.data = rankedData;
	}
	
	/**
	 * Sequential step 2:
	 * @param fasta
	 * @param decoyPrefix
	 */
	private void assignProteinId (FastaLoader fasta, String decoyPrefix) {
		Hashtable<String, Boolean> check = new Hashtable<String, Boolean>();
		ArrayList<String> sequences = new ArrayList<String>();
		for(HarmonyData hData : this.data) {
			if(check.get(hData.sequence) == null) {
				sequences.add(hData.sequence);
				check.put(hData.sequence, true);
			}
		}
		
		Pattern geneNamePattern = Pattern.compile("(GN=[\\S]*)");
		Pattern proteinEvidenceGroup = Pattern.compile("(PE=[0-9]*)");
		
		Trie trie = Trie.builder().addKeywords(sequences).build();
		Hashtable<String, ArrayList<HarmonyProtein>> seqToProteins = new Hashtable<String, ArrayList<HarmonyProtein>>();
		for(FastaEntry entry : fasta.entries) {
			String header = entry.originHeader;
			String sequence = entry.sequence;
			
			// skip decoy
			if(header.startsWith(decoyPrefix)) {
				continue;
			}
			
			Collection<Emit> emits = trie.parseText(sequence);
			
			if(emits.size() > 0) {
				// find GN (gene name)
				Matcher matcher = geneNamePattern.matcher(header);
				String geneName = "Unknown";
				int proteinEvidence = 999;
				if(matcher.find()) {
					geneName = matcher.group().replace("GN=", "");
				}
				
				matcher = proteinEvidenceGroup.matcher(header);
				if(matcher.find()) {
					proteinEvidence = Integer.parseInt(matcher.group().replace("PE=", ""));
				}
				
				
				String[] split = header.split("\\s")[0].split("\\|");
				header = split[0];
				for(int i=1; i<split.length; i++) {
					header += ":"+split[i];
				}
				
				for(Emit emit : emits) {
					ArrayList<HarmonyProtein> proteins = seqToProteins.get(emit.getKeyword());
					if(proteins == null) {
						proteins = new ArrayList<HarmonyProtein>();
					}
					
					proteins.add(new HarmonyProtein(header, geneName, proteinEvidence));
					seqToProteins.put(emit.getKeyword(), proteins);
				}
			}
		}
		
		// select better protein evidence (lower is better)
		Hashtable<String, Integer> minPEs = new Hashtable<String, Integer>();
		seqToProteins.forEach((sequence, proteins)->{
			Collections.sort(proteins);
			
			// find minimum PE
			int minPE = 999;
			for(HarmonyProtein hProtein : proteins) {
				minPE = Math.min(hProtein.proteinEvidence, minPE);
			}
			// save min PE
			minPEs.put(sequence, minPE);
			
			// remain only minimum PE
			int size = proteins.size();
			for(int i=size-1; i>=0 ;i--) {
				if(proteins.get(i).proteinEvidence != minPE) {
					proteins.remove(i);
				}
			}
		});
		
		
		Hashtable<String, String> uniqueProteins = new Hashtable<String, String>();
		Hashtable<String, String> uniqueGeneNames = new Hashtable<String, String>();
		// remove duplications
		seqToProteins.forEach((sequence, proteins)->{
			check.clear();
			// prevent to add Unknown
			check.put("Unknown", true); 
			StringBuilder pIds = new StringBuilder();
			StringBuilder gIds = new StringBuilder();
			for(HarmonyProtein hProtein : proteins) {
				if(check.get(hProtein.proteinId) == null) {
					pIds.append(hProtein.proteinId).append("|");
					check.put(hProtein.proteinId, true);
				}
				
				if(check.get(hProtein.geneName) == null) {
					gIds.append(hProtein.geneName).append("|");
					check.put(hProtein.geneName, true);
				}
			}
			
			// delete at the end of "|"
			if(pIds.length() == 0) {
				pIds.append("Unknown");
			} else {
				pIds.deleteCharAt(pIds.length()-1);
			}
			
			if(gIds.length() == 0) {
				gIds.append("Unknown");
			} else {
				gIds.deleteCharAt(gIds.length()-1);
			}
						
			uniqueProteins.put(sequence, pIds.toString());
			uniqueGeneNames.put(sequence, gIds.toString());
		});
		
		// add proteins and genes
		for(HarmonyData hData : this.data) {
			String proteinId = uniqueProteins.get(hData.sequence);
			String geneName = uniqueGeneNames.get(hData.sequence);
			Integer minPE = minPEs.get(hData.sequence);
			
			if(proteinId == null) {
				proteinId = "Unknown"; // empty
			}
			
			if(geneName == null) {
				geneName = "Unknown"; // empty
			}
			
			if(minPE == null) {
				minPE = 999;
			}
			
			hData.proteinId = proteinId;
			hData.proteinGeneName = geneName;
			hData.proteinEvidence = minPE;
			hData.setMajorAnnotation();
		}
		
	}
	

	/**
	 * Sequential step 3:
	 * Calculate FDR for high (1%) and medium (5%) confidence level.
	 * 
	 * 
	 */
	private void assignConfidence () {
		Hashtable<String, Boolean> classes = new Hashtable<String, Boolean>();
		// find classes
		for(HarmonyData hData : this.data) {
			classes.put(hData.fdrClassId, true);
		}
		// check classes
		System.out.println("Detected classes: "+classes.size());
		// sort by score
		Collections.sort(this.data, new Comparator<HarmonyData>() {
			@Override
			public int compare(HarmonyData o1, HarmonyData o2) {
				double s1 = Double.parseDouble(o1.score);
				double s2 = Double.parseDouble(o2.score);
				if(s1 < s2) {
					return 1;
				} else if(s1 > s2) {
					return -1;
				}
				return 0;
			}
			
		});
		// FDR
		classes.forEach((class_, nil)->{
			double fdr = 0;
			double t = 0;
			double d = 0;
			int fdr001Idx = -1;
			int fdr005Idx = -1;
			
			for(int i=0; i< this.data.size(); i++) {
				HarmonyData hData  = this.data.get(i);
				if(hData.fdrClassId.equals(class_)) {
					if(hData.label.equalsIgnoreCase("1")) {
						t++;
						
						fdr = d/t;
						
						if(fdr < 0.01) {
							fdr001Idx = i;
						}
						
						if(fdr < 0.05) {
							fdr005Idx = i;
						}
						
					} else {
						d++;
					}
					
				}
			}
			
			// assign
			for(int i=0; i< this.data.size(); i++) {
				HarmonyData hData  = this.data.get(i);
				if(hData.fdrClassId.equals(class_)) {
					// fdr 0.01
					if(i <= fdr001Idx) {
						hData.confidence = "High";
					} else if(i <= fdr005Idx) {
						hData.confidence = "Medium";
					} 
				}
			}
		});
		
		// remove null confidence
		ArrayList<HarmonyData> confidentData = new ArrayList<HarmonyData>();
		for(int i=0; i< this.data.size(); i++) {
			HarmonyData hData  = this.data.get(i);
			if(hData.confidence != null) {
				confidentData.add(hData);
			}
		}
		
		this.data = confidentData;
	}

	/**
	 * Sequential step 4:
	 */
	private void dropDecoyLabel () {
		ArrayList<HarmonyData> tData = new ArrayList<HarmonyData>();
		for(int i=0; i<this.data.size(); i++) {
			if(this.data.get(i).label.equalsIgnoreCase("1")) {
				tData.add(this.data.get(i));
			}
		}
		this.data = tData;
	}
	
	public void write(File file) throws IOException {
		BufferedWriter BW = new BufferedWriter(new FileWriter(file));

		for(int i=0; i<HEADER.length; i++) {
			if(i!=0) {
				BW.append("\t");
			}
			BW.append(HEADER[i]);
		}
		BW.newLine();
		for(int i=0; i<data.size(); i++) {
			BW.append(data.get(i).toString());
			BW.newLine();
		}
		
		BW.close();
	}
	
	
	
	/**
	 * Case-insensitive.
	 * 
	 * @param taget
	 * @param replacement
	 */
	public void renameHeader (String target, String replacement) {
		for(int i=0; i<HEADER.length; i++) {
			if(HEADER[i].equalsIgnoreCase(target)) {
				HEADER[i] = replacement;
			}
		}
	}
	
	public static HarmonyResult merge (HarmonyResult[] results) {
		HarmonyResult mergedResult = results[0];
		
		Hashtable<String, ArrayList<HarmonyData>> aggregateByPSMKey = new Hashtable<String, ArrayList<HarmonyData>>();
		
		// aggregate by PSM
		for(int i=0; i<results.length; i++) {
			HarmonyResult result = results[i];
			for(HarmonyData hData : result.data) {
				ArrayList<HarmonyData> hList = aggregateByPSMKey.get(hData.psmKey);
				if(hList == null) {
					hList = new ArrayList<HarmonyData>();
				}
				hList.add(hData);
				aggregateByPSMKey.put(hData.psmKey, hList);
			}
		}
		
		// initialize
		mergedResult.data.clear();
		
		// for each PSM, merge their RNA reads
		aggregateByPSMKey.forEach((psmKey, hList)->{
			Hashtable<String, ArrayList<HarmonyData>> aggregateByGenomicKey = new Hashtable<String, ArrayList<HarmonyData>>();
			for(HarmonyData hData : hList) {
				ArrayList<HarmonyData> gList = aggregateByGenomicKey.get(hData.genomicKey);
				if(gList == null) {
					gList = new ArrayList<HarmonyData>();
				}
				gList.add(hData);
				aggregateByGenomicKey.put(hData.genomicKey, gList);
				
				
			}
			
			// sum of reads
			// geometric mean of RPHM and then recalculate proportion.
			// 
			aggregateByGenomicKey.forEach((genomicKey, gList)->{
				int sumReads = 0;
				double avgRPHM = 0;
				
				for(HarmonyData gData : gList) {
					sumReads += Integer.parseInt(gData.matchedReadCount);
					
					double rphm = Double.parseDouble(gData.matchedRPHM);
					if(rphm > 0) {
						if(avgRPHM == 0) {
							avgRPHM = rphm;
						} else {
							avgRPHM *= rphm;
						}
					}
				}
				
				avgRPHM = Math.pow(avgRPHM, 1.0/results.length);
				
				// select random harmony data. Indeed, they are exactly same except for RNA abundance.
				HarmonyData gData = gList.get(0);
				// update abundance
				gData.matchedReadCount = sumReads+"";
				gData.matchedRPHM = avgRPHM+"";
				mergedResult.data.add(gData);
			});
			
		});
		

		// aggregate by PSM
		aggregateByPSMKey.clear();
		for(HarmonyData hData : mergedResult.data) {
			ArrayList<HarmonyData> hList = aggregateByPSMKey.get(hData.psmKey);
			if(hList == null) {
				hList = new ArrayList<HarmonyData>();
			}
			hList.add(hData);
			aggregateByPSMKey.put(hData.psmKey, hList);
		}
		
		// recalculate proportion by RPHM values
		aggregateByPSMKey.forEach((psmKey, hList)->{
			double sumOfRPHM = 0;
			for(HarmonyData hData : hList) {
				sumOfRPHM += Double.parseDouble(hData.matchedRPHM);
			}
			
			for(HarmonyData hData : hList) {
				if(sumOfRPHM == 0) {
					hData.proportion = "0";
				} else {
					hData.proportion = (Double.parseDouble(hData.matchedRPHM) / sumOfRPHM)+"";
				}
			}
		});
		

		System.out.println(mergedResult.data.size());
		
		return mergedResult;
	}
}
