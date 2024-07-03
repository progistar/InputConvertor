package zhanglab.inputconvertor.module;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.data.FastaEntry;
import zhanglab.inputconvertor.data.FastaLoader;
import zhanglab.inputconvertor.data.Peptide;
import zhanglab.inputconvertor.env.InputConvertorConstants;

class CometXML {
	public String title = null;
	public String peptide = null;
}

class CometTXT {
	public String originRecord = null;
	public String specId = null;
	public String genomicId = null;
	public String label = null;
	public String icTitle = null;
	public String icScanNum = null;
	public String icObRT = null;
	public String icCharge = null;
	public String icScore = null;
	public String icPeptide = null;
	public String inferredPeptide = null;
	public String proteinId = null;
	public String isCanonical = "true";

	public String toString() {
		StringBuilder str = new StringBuilder();
		
		str.append(specId).append("\t")
		.append(label).append("\t")
		.append(icTitle).append("\t")
		.append(icScanNum).append("\t")
		.append(icObRT).append("\t")
		.append(icCharge).append("\t")
		.append(icScore).append("\t")
		.append(icPeptide).append("\t")
		.append(inferredPeptide).append("\t")
		.append(genomicId).append("\t")
		.append(isCanonical).append("\t")
		.append(originRecord);
		
		return str.toString();
	}
}

public class ConvertCometTopXg {

	public static String inputTXTFilePath;
	public static String inputPINFilePath;
	public static String inputXMLFilePath;
	public static String outputFilePath;
	public static String fastaFilePath;
	public static Hashtable<String, String> targetClassHash = new Hashtable<String, String>();
	public static Hashtable<String, String> decoyClassHash = new Hashtable<String, String>();
	public static String decoyPrefix = "XXX";
	
	public static void convertTopXgOutput (String[] args) throws IOException {
		parseOptions(args);
		
		int modificationIdx = -1;
		int proteinIdx = -1;
		int scanIdx = -1;
		int chargeIdx = -1;
		int xcorrIdx = -1;
		int rtIdx = -1;
		File inputTXTFile = new File(inputTXTFilePath);
		File inputPINFile = new File(inputPINFilePath);
		File inputXMLFile = new File(inputXMLFilePath);
		File outputFile = new File(outputFilePath);
		
		ArrayList<CometTXT> cometTXTs = new ArrayList<CometTXT>();
		ArrayList<CometXML> cometXMLs = loadCometXML(inputXMLFile);
		
		boolean isAppended = outputFile.exists();
		
		BufferedReader BR = new BufferedReader(new FileReader(inputTXTFile));
		BufferedReader PIN= new BufferedReader(new FileReader(inputPINFile));
		BufferedWriter BW = new BufferedWriter(new FileWriter(outputFile, isAppended));
		
		String line = null;
		
		BR.readLine(); //skip meta
		////////////////////////////////// Build header ////////////////////////////////////////////////
		String header = BR.readLine();
		String[] headerSplit = header.split("\t");
		modificationIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.COMET_MODIFICATION_FIELD_NAME);
		proteinIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.COMET_PROTEIN_FIELD_NAME);
		scanIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.COMET_SCAN_FIELD_NAME);
		chargeIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.COMET_CHARGE_FIELD_NAME);
		xcorrIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.COMET_XCORR_FIELD_NAME);
		rtIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.COMET_RT_FIELD_NAME);
		////////////////////////////// End of header ////////////////////////////////////////////////
		
		/////// PIN READER /////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////Build header ////////////////////////////////////////////////
		String pinHeader = PIN.readLine();
		headerSplit = pinHeader.split("\t");
		int pinProteinIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.COMET_PIN_PROTEIN_FIELD_NAME);
		int pinPeptideIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.COMET_PIN_PEPTIDE_FIELD_NAME);
		int pinSpecIdIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.COMET_PIN_SPECID_FIELD_NAME);
		int pinPeplenIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.COMET_PIN_PEPLEN_FIELD_NAME);
		////////////////////////////// End of header ////////////////////////////////////////////////
		
		
		if(!isAppended) {
			header = InputConvertorConstants.IC_SPEC_ID_FEILD_NAME +"\t"
					+ InputConvertorConstants.IC_LABEL_FEILD_NAME +"\t"
					+ InputConvertorConstants.IC_TITLE_FIELD_NAME +"\t"
					+ InputConvertorConstants.IC_SCAN_NUM_FIELD_NAME +"\t"
					+ InputConvertorConstants.IC_RT_FIELD_NAME +"\t"
					+ InputConvertorConstants.IC_CHARGE_FIELD_NAME +"\t"
					+ InputConvertorConstants.IC_SEARCH_SCORE_FIELD_NAME +"\t"
					+ InputConvertorConstants.IC_PEPTIDE_FIELD_NAME +"\t"
					+ InputConvertorConstants.IC_INFERRED_PEPTIDE_FIELD_NAME +"\t"
					
					+ InputConvertorConstants.IC_GENOMIC_ID_FEILD_NAME +"\t"
					+ InputConvertorConstants.IC_IS_CANONICAL_FEILD_NAME +"\t"
					
					+ header;
			
			
			header += "\t###\t" + pinHeader.replace("PepLen", "Length=7\tLength=8\tLength>8");
			BW.append(header);
			BW.newLine();
		}
		
		Hashtable<String, Integer> genomicIdx = new Hashtable<String, Integer>();
		int idx = 0;

		while((line = BR.readLine()) != null) {
			
			String[] fields = line.split("\t");
			CometTXT cometTXT = new CometTXT();
			CometXML cometXML = cometXMLs.get(idx++);
			
			cometTXT.icTitle = cometXML.title;
			cometTXT.icPeptide = cometXML.peptide;
			cometTXT.icScanNum = fields[scanIdx];
			cometTXT.icCharge = fields[chargeIdx];
			cometTXT.icScore = fields[xcorrIdx];
			cometTXT.proteinId = fields[proteinIdx];
			cometTXT.icObRT = fields[rtIdx];
			cometTXT.originRecord = line;
			
			String[] modifications = fields[modificationIdx].split("\\,");
			String[] modLoc = new String[cometTXT.icPeptide.length()+1];
			for(int i=0; i<modifications.length; i++) {
				String modification = modifications[i];
				if(modification.equalsIgnoreCase("-")) continue;
				String[] modInfo = modification.split("\\_");
				Integer loc = Integer.parseInt(modInfo[0]);
				String mod = modInfo[2];
				modLoc[loc] = mod;
			}
			
			StringBuilder modPeptide = new StringBuilder();
			// N-term mod
			if(modLoc[0] != null) {
				modPeptide.append("[").append(modLoc[0]).append("]");
			}
			// Mid and C-term mods
			for(int i=0; i<cometTXT.icPeptide.length(); i++) {
				modPeptide.append(cometTXT.icPeptide.charAt(i));
				if(modLoc[i+1] != null) {
					modPeptide.append("[").append(modLoc[i+1]).append("]");
				}
			}
			Peptide peptide = new Peptide(modPeptide.toString(), InputConvertorConstants.COMET);
			cometTXT.icPeptide = peptide.modPeptide;
			cometTXT.inferredPeptide = peptide.stripPeptide;
			
			// label setting
			if(cometTXT.proteinId.startsWith(decoyPrefix)) {
				cometTXT.label = "-1";
				cometTXT.isCanonical = decoyClassHash.get(cometXML.peptide.replace("I", "L"));
			} else {
				cometTXT.label = "1";
				cometTXT.isCanonical = targetClassHash.get(cometXML.peptide.replace("I", "L"));
			}
			if(cometTXT.isCanonical == null) {
				cometTXT.isCanonical = "true";
			}
			
			// sepcId
			cometTXT.specId = cometTXT.icTitle +"|" +cometTXT.icScanNum +"|"+cometTXT.icCharge;
			
			cometTXTs.add(cometTXT);
			
			// corresponding PIN file
			String pinRecord = PIN.readLine();
			String[] pinFields = pinRecord.split("\t");
			StringBuilder pinStr = new StringBuilder();
			
			pinFields[pinSpecIdIdx] = cometTXT.specId;
			pinFields[pinPeptideIdx] = cometTXT.icPeptide;
			int pepLen = Integer.parseInt(pinFields[pinPeplenIdx]);
			
			String pepLenFeature = "0\t0\t1";
			if(pepLen == 7) {
				pepLenFeature = "1\t0\t0";
			} else if(pepLen == 8) {
				pepLenFeature = "0\t1\t0";
			}
			pinFields[pinPeplenIdx] = pepLenFeature;
			
			// genomic Id
			String key = pinFields[pinPeptideIdx]+"_"+pinFields[pinProteinIdx];
			Integer gIdx = genomicIdx.get(key);
			if(gIdx == null) {
				gIdx = genomicIdx.size()+1;
				genomicIdx.put(key, gIdx);
			}
			
			if(cometTXT.label.equalsIgnoreCase("-1")) {
				pinFields[pinProteinIdx] = decoyPrefix+"_"+gIdx;
			} else{
				pinFields[pinProteinIdx] = ""+gIdx;
			}
			cometTXT.genomicId = ""+gIdx;
			
			pinStr.append(cometTXT.originRecord).append("\t###");
			for(int i=0; i<=pinProteinIdx; i++) {
				pinStr.append("\t").append(pinFields[i]);
			}
			
			cometTXT.originRecord = pinStr.toString();
			
			// skip selenocysteine!
			if(cometTXT.icPeptide.contains("U")) {
				continue;
			}
			
			BW.append(cometTXT.toString());
			BW.newLine();
		}
		
		BR.close();
		BW.close();
		PIN.close();
		
	}
	
	public static ArrayList<CometXML> loadCometXML (File file) throws IOException {
		ArrayList<CometXML> xmls = new ArrayList<CometXML>();
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		CometXML xml = null;
		String title = null;
		
		Hashtable<String, Boolean> isIncluded = new Hashtable<String, Boolean>();
		ArrayList<String> peptides = new ArrayList<String>();
		while((line = BR.readLine()) != null) {
			line = line.trim();
			String[] fields = null;
			if(line.startsWith("<spectrum_query spectrum=")) {
				// spectrumNatvieID is inconsistent between mzML and mgf.
				// spectrum is consistent; however, its scan number is not zero padded.
				fields = line.split("spectrum=");
				String tmpTitle = fields[1].split("\\s")[0].replace("\"", "");
				// correct scanNum
				String[] tmpTitleArray = tmpTitle.split("\\.");
				String scanNum = tmpTitleArray[tmpTitleArray.length-2];
				int scanInt = Integer.parseInt(scanNum);
				tmpTitleArray[tmpTitleArray.length-2] = scanInt+"";
				tmpTitleArray[tmpTitleArray.length-3] = scanInt+"";
				
				title = "";
				for(int i=0; i<tmpTitleArray.length; i++) {
					if(i!=0) {
						title += ".";
					}
					title += tmpTitleArray[i];
				}
			} else if(line.startsWith("<search_hit hit_rank=")) {
				
				fields = line.split("\\s");
				for(int i=0; i<fields.length; i++) {
					if(fields[i].startsWith("peptide=")) {
						String peptide = fields[i].split("\\=")[1].trim().replace("\"", "");
						xml = new CometXML();
						xmls.add(xml);
						xml.title = title;
						xml.peptide = peptide;
						
						// I/L replace
						peptide = peptide.replace("I", "L");
						if(isIncluded.get(peptide) == null) {
							isIncluded.put(peptide, true);
							peptides.add(peptide);
						}
					}
				}
			}
		}
		
		BR.close();
		

		// load fasta file
		FastaLoader fastaLoader = new FastaLoader(new File(fastaFilePath));
		Trie trie = Trie.builder().addKeywords(peptides).build();
		
		ArrayList<FastaEntry> entries = fastaLoader.entries;
		for(FastaEntry entry : entries) {
			StringBuilder proteinSequence = new StringBuilder(entry.sequence.replace("I", "L"));
			boolean isTarget = true;
			if(entry.originHeader.startsWith(decoyPrefix)) {
				isTarget = false;
			}
			
			Collection<Emit> emits = trie.parseText(proteinSequence);
			if(emits.size() == 0) continue;
			
			String class_ = "true";
			if(entry.originHeader.contains(InputConvertorConstants.NON_REF_HEADER_ID)) {
				class_ = "false";
			}
			
			for(Emit emit : emits) {
				String peptide = emit.getKeyword();
				String thatClass = null;
				String thisClass = "true";
				if(isTarget) {
					thatClass = targetClassHash.get(peptide);
				} else {
					thatClass = decoyClassHash.get(peptide);
				}
				
				if(thatClass == null) {
					thatClass = class_;
				}
				
				if(thatClass.equalsIgnoreCase("true") || class_.equalsIgnoreCase("true")) {
					thisClass = "true";
				} else {
					thisClass = "false";
				}
				
				if(isTarget) {
					targetClassHash.put(peptide, thisClass);
				}else {
					decoyClassHash.put(peptide, thisClass);
				}
			}
		}
		
		
		return xmls;
	}

	public static void parseOptions (String[] args) {
		CommandLine cmd = null;
		Options options = new Options();
		
		// Mandatory
		Option optionInputTxt = Option.builder("i")
				.longOpt("input_txt").argName("txt")
				.hasArg()
				.required(true)
				.desc("text output of comet result")
				.build();
		
		Option optionInputPin = Option.builder("p")
				.longOpt("input_pin").argName("pin")
				.hasArg()
				.required(true)
				.desc("pin output of comet result")
				.build();
		
		Option optionInputXML = Option.builder("x")
				.longOpt("input_xml").argName("xml")
				.hasArg()
				.required(true)
				.desc("XML output of comet result")
				.build();
		
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("tsv")
				.hasArg()
				.required(true)
				.desc("output")
				.build();
		
		Option optionFasta = Option.builder("f")
				.longOpt("fasta").argName("fa|fasta")
				.hasArg()
				.required(true)
				.desc("sequence database")
				.build();
		
		options.addOption(optionInputTxt)
		.addOption(optionInputPin)
		.addOption(optionInputXML)
		.addOption(optionFasta)
		.addOption(optionOutput);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
	    ArrayList<String> tmpArgs = new ArrayList<String>();
	    for(int i=0; i<args.length; i++) {
	    	if(args[i].equalsIgnoreCase("-i") || args[i].equalsIgnoreCase("--input_txt") ||
			args[i].equalsIgnoreCase("-o") || args[i].equalsIgnoreCase("--output") ||
			args[i].equalsIgnoreCase("-p") || args[i].equalsIgnoreCase("--input_pin") ||
			args[i].equalsIgnoreCase("-x") || args[i].equalsIgnoreCase("--input_xml") ||
			args[i].equalsIgnoreCase("-f") || args[i].equalsIgnoreCase("--fasta")) {
	    		tmpArgs.add(args[i++]);
	    		tmpArgs.add(args[i]);
	    	}
	    }
	    
	    String[] nArgs = new String[tmpArgs.size()];
	    for(int i =0; i<tmpArgs.size(); i++) {
	    	nArgs[i] = tmpArgs.get(i);
	    }
	    
	    
		try {
		    cmd = parser.parse(options, nArgs, false);
		    
		    if(cmd.hasOption("i")) {
		    	inputTXTFilePath = cmd.getOptionValue("i");
		    }
		    
		    if(cmd.hasOption("p")) {
		    	inputPINFilePath = cmd.getOptionValue("p");
		    }
		    
		    if(cmd.hasOption("x")) {
		    	inputXMLFilePath = cmd.getOptionValue("x");
		    }
		    
		    if(cmd.hasOption("o")) {
		    	outputFilePath = cmd.getOptionValue("o");
		    }
		    
		    if(cmd.hasOption("f")) {
		    	fastaFilePath = cmd.getOptionValue("f");
		    }
		    
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			isFail = true;
		}
		
		if(isFail) {
		    helper.printHelp("Usage:", options);
		    System.exit(0);
		}
		
		System.out.println();
	}
}
