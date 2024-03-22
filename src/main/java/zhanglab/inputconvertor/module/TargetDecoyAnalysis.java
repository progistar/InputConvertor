package zhanglab.inputconvertor.module;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.data.Peptide;
import zhanglab.inputconvertor.env.InputConvertorConstants;

class TDRecord implements Comparable<TDRecord>{
	String record;
	double score;
	int label;
	int peptideLength;
	double fdr;
	
	@Override
	public int compareTo(TDRecord o) {
		if(this.score < o.score) {
			return 1;
		} else if(this.score > o.score) {
			return -1;
		} else if(this.label > o.label) {
			return -1;
		} else if(this.label < o.label) {
			return 1;
		}
		return 0;
	}
	
	public String toString() {
		return record + "\t" + score;
	}
}

public class TargetDecoyAnalysis {

	public static String inputFilePath = null;
	public static String outputFilePath = null;
	public static String targetFilePath = null;
	public static String decoyFilePath = null;
	public static double fdr = -1;
	public static double minScore = 0;
	
	private TargetDecoyAnalysis() {};
	
	
	public static void doFDR (String[] args) throws IOException, ParseException {
		
        parseOptions(args);
        
        Hashtable<String, String[]> targetPSMs = getHashTableFromPSM(targetFilePath);
        Hashtable<String, String[]> decoyPSMs = getHashTableFromPSM(decoyFilePath);
        
        
        BufferedReader BR = new BufferedReader(new FileReader(inputFilePath));
        BufferedWriter BW = new BufferedWriter(new FileWriter(outputFilePath));
        String line = null;
        String header = BR.readLine();
        
        // build header ////////////////////////////////////////////////
        BW.append(header).append("\t")
        .append("percolator_score");
        BW.newLine();
        ////////////////////////////////////// End of building header /
        
        String[] headerSplit = header.split("\t");
        int specIdIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.IC_SPEC_ID_FEILD_NAME);
        int labelIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.IC_LABEL_FEILD_NAME);
        int genomicIdIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.IC_GENOMIC_ID_FEILD_NAME);
        int isCanonicalIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.IC_IS_CANONICAL_FEILD_NAME);
        int icPeptideIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.IC_PEPTIDE_FIELD_NAME);
        int inferredPeptideIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.IC_INFERRED_PEPTIDE_FIELD_NAME);

        
        ArrayList<TDRecord> cRecords = new ArrayList<TDRecord>();
        ArrayList<TDRecord> ncRecords = new ArrayList<TDRecord>();
        
        Hashtable<String, String> duplications = new Hashtable<String, String>();
        while((line = BR.readLine()) != null) {
        	String[] fields = line.split("\t");
        	String specId = fields[specIdIdx];
        	String genomicId = fields[genomicIdIdx];
        	String isCanonical = fields[isCanonicalIdx];
        	String label = fields[labelIdx];
        	Peptide peptide = new Peptide(fields[icPeptideIdx], InputConvertorConstants.IC_CONVERTOR);
        	peptide.icPeptideToILDetermination(fields[inferredPeptideIdx]);
        	int peptideLength = fields[inferredPeptideIdx].length();
        	
        	// key: specId + peptide.modPeptide
        	// this is because same genomic id and spec id can be occurred when there are modified peptides (site-localization problem).
        	// Casanovo sometimes distinguishes I/L... 
        	String key = specId+"_"+peptide.modPeptide;
        	
        	String[] psm = targetPSMs.get(key) != null ? targetPSMs.get(key) : decoyPSMs.get(key);
        	if(psm == null) {
        		continue;
        	}
        	
        	String genomicIdInPSM = psm[0].replace("XXX_", "");
        	
        	if(genomicIdInPSM.equalsIgnoreCase(genomicId)) {
        		TDRecord record = new TDRecord();
        		record.record = line;
        		record.score = Double.parseDouble(psm[1]);
        		record.label = Integer.parseInt(label);
        		record.peptideLength = peptideLength;
        		
        		if(isCanonical.equalsIgnoreCase("true")) {
        			cRecords.add(record);
        		} else {
        			ncRecords.add(record);
        		}
        	}
        }
        
        // length-specific FDR control
        ArrayList<TDRecord> passTargetList = new ArrayList<TDRecord>();
        boolean isLengthSpecific = false;
        if(isLengthSpecific) {
        	for(int len=7; len<=9; len++) {
            	int startLen = len;
            	int endLen = len == 9 ? 100 : len;
            	
            	ArrayList<TDRecord> lengthSpecificCRecords = new ArrayList<TDRecord>();
            	ArrayList<TDRecord> lengthSpecificNCRecords = new ArrayList<TDRecord>();
            	
            	
            	for(TDRecord record : cRecords) {
            		if(startLen <= record.peptideLength && record.peptideLength <= endLen) {
            			lengthSpecificCRecords.add(record);
            		}
            	}
            	
            	for(TDRecord record : ncRecords) {
            		if(startLen <= record.peptideLength && record.peptideLength <= endLen) {
            			lengthSpecificNCRecords.add(record);
            		}
            	}
            	
            	System.out.println("Canonical PSMs with length "+len);
            	System.out.println("Target PSMs: "+getNumOfRecords(lengthSpecificCRecords, 1));
            	System.out.println("Decoy PSMs: "+getNumOfRecords(lengthSpecificCRecords, -1));
            	passTargetList.addAll(getFDR(lengthSpecificCRecords, fdr));
            	System.out.println("Non-canonical PSMs with length "+len);
            	System.out.println("Target PSMs: "+getNumOfRecords(lengthSpecificNCRecords, 1));
            	System.out.println("Decoy PSMs: "+getNumOfRecords(lengthSpecificNCRecords, -1));
            	passTargetList.addAll(getFDR(lengthSpecificNCRecords, fdr));
            }
        } else {
        	int len7=0; int len8=0; int len9=0;
        	System.out.println("Canonical PSMs");
        	System.out.println("Target PSMs: "+getNumOfRecords(cRecords, 1));
        	System.out.println("Decoy PSMs: "+getNumOfRecords(cRecords, -1));
        	ArrayList<TDRecord> pass = getFDR(cRecords, fdr);
        	for(int i=0; i<pass.size(); i++) {
        		if(pass.get(i).label > 0 && pass.get(i).peptideLength == 7) {
        			len7++;
        		} else if(pass.get(i).label > 0 && pass.get(i).peptideLength == 8) {
        			len8++;
        		} else if(pass.get(i).label > 0 && pass.get(i).peptideLength >= 9) {
        			len9++;
        		}
        	}
        	System.out.println("IDs: "+len7+"|"+len8+"|"+len9);
        	passTargetList.addAll(pass);
        	System.out.println("Non-canonical PSMs");
        	System.out.println("Target PSMs: "+getNumOfRecords(ncRecords, 1));
        	System.out.println("Decoy PSMs: "+getNumOfRecords(ncRecords, -1));
        	pass = getFDR(ncRecords, fdr);
        	len7=0; len8=0; len9=0;
        	for(int i=0; i<pass.size(); i++) {
        		if(pass.get(i).label > 0 && pass.get(i).peptideLength == 7) {
        			len7++;
        		} else if(pass.get(i).label > 0 && pass.get(i).peptideLength == 8) {
        			len8++;
        		} else if(pass.get(i).label > 0 && pass.get(i).peptideLength >= 9) {
        			len9++;
        		}
        	}
        	System.out.println("IDs: "+len7+"|"+len8+"|"+len9);
        	passTargetList.addAll(pass);
        }
        
        
        for(int i=0; i<passTargetList.size(); i++) {
        	TDRecord record = passTargetList.get(i);
        	
        	String specId = record.record.split("\t")[specIdIdx];
        	if(duplications.get(specId) == null) {
        		duplications.put(specId, specId);
        		BW.append(record.toString());
            	BW.newLine();
        	} else {
        		System.out.println("Duplicated records: "+specId);
        	}
        }
        
        BW.close();
        BR.close();
	}
	
	private static int getNumOfRecords (ArrayList<TDRecord> records, int label) {
		int tCount = 0;
		int dCount = 0;
		for(int i=0; i<records.size(); i++) {
			TDRecord record = records.get(i);
			
			if(record.label > 0) {
				tCount++;
			} else {
				dCount++;
			}
		}
		if(label > 0) {
			return tCount;
		} else {
			return dCount;
		}
	}
	
	private static ArrayList<TDRecord> getFDR (ArrayList<TDRecord> records, double fdr) {
		Collections.sort(records);
		
		double tCount = 0;
		double dCount = 0;
		double estimatedFDR = 0;
		int lastIdx   = -1;
		for(int i=0; i<records.size(); i++) {
			TDRecord record = records.get(i);
			
			if(record.label > 0) {
				tCount++;
			} else {
				dCount++;
			}
			
			if(tCount == 0) {
				record.fdr = 1;
			} else {
				record.fdr = dCount/tCount;
				
				if(record.fdr < fdr) {
					lastIdx = i;
					estimatedFDR = record.fdr;
				}
			}
		}
		
		ArrayList<TDRecord> passTargetList = new ArrayList<TDRecord>();
		for(int i=0; i<=lastIdx; i++) {
			TDRecord record = records.get(i);
			if(record.label > 0) {
				passTargetList.add(record);
			}
		}
		
		System.out.println("FDR: " + estimatedFDR);
		System.out.println("IDs: "+passTargetList.size());
		
		return passTargetList;
	}
	
	
	
	private static Hashtable<String, String[]> getHashTableFromPSM (String fileName) throws IOException {
		File file = new File(fileName);
		Hashtable<String, String[]> table = new Hashtable<String, String[]>();
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		String[] header = BR.readLine().split("\t");
		int psmIdIdx = InputConvertorConstants.getFieldIndex(header, "PSMId");
		int scoreIdx = InputConvertorConstants.getFieldIndex(header, "score");
		int proteinIdIdx = InputConvertorConstants.getFieldIndex(header, "proteinIds");
		int peptideIdx = InputConvertorConstants.getFieldIndex(header, "peptide");
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String[] record = {fields[proteinIdIdx], fields[scoreIdx]};
			String key = fields[psmIdIdx]+"_"+fields[peptideIdx];
			table.put(key, record);
		}
		
		BR.close();
		
		return table;
	}
	
	public static void parseOptions (String[] args) {
		CommandLine cmd = null;
		Options options = new Options();
		
		// Mandatory
		Option optionInput = Option.builder("i")
				.longOpt("input").argName("pXg")
				.hasArg()
				.required(true)
				.desc("pXg result")
				.build();
		
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("tsv")
				.hasArg()
				.required(true)
				.desc("fdr output file")
				.build();
		
		Option optionTargetFile = Option.builder("t")
				.longOpt("target").argName("psm")
				.hasArg()
				.required(true)
				.desc("target psm file from Percolator")
				.build();
		
		Option optionDecoyFile = Option.builder("d")
				.longOpt("decoy").argName("psm")
				.hasArg()
				.required(true)
				.desc("decoy psm file from Percolator")
				.build();
		
		Option optionFDR = Option.builder("f")
				.longOpt("fdr").argName("float")
				.hasArg()
				.required(true)
				.desc("FDR value")
				.build();
		
		Option optionMinScore = Option.builder("s")
				.longOpt("score").argName("float")
				.hasArg()
				.required(false)
				.desc("Min percolator score (default: 0)")
				.build();
		
		options.addOption(optionInput)
		.addOption(optionTargetFile)
		.addOption(optionDecoyFile)
		.addOption(optionFDR)
		.addOption(optionOutput)
		.addOption(optionMinScore);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
	    ArrayList<String> tmpArgs = new ArrayList<String>();
	    for(int i=0; i<args.length; i++) {
	    	if(args[i].equalsIgnoreCase("-i") || args[i].equalsIgnoreCase("--input") ||
			args[i].equalsIgnoreCase("-t") || args[i].equalsIgnoreCase("--target") ||
			args[i].equalsIgnoreCase("-d") || args[i].equalsIgnoreCase("--decoy") ||
			args[i].equalsIgnoreCase("-f") || args[i].equalsIgnoreCase("--fdr") ||
			args[i].equalsIgnoreCase("-o") || args[i].equalsIgnoreCase("--output") ||
			args[i].equalsIgnoreCase("-s") || args[i].equalsIgnoreCase("--score")) {
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
		    	inputFilePath = cmd.getOptionValue("i");
		    }
		    
		    if(cmd.hasOption("o")) {
		    	outputFilePath = cmd.getOptionValue("o");
		    }
		    
		    if(cmd.hasOption("d")) {
		    	decoyFilePath = cmd.getOptionValue("d");
		    }
		    
		    if(cmd.hasOption("t")) {
		    	targetFilePath = cmd.getOptionValue("t");
		    }
		    
		    if(cmd.hasOption("f")) {
		    	fdr = Double.parseDouble(cmd.getOptionValue("f"));
		    }
		    
		    if(cmd.hasOption("s")) {
		    	minScore = Double.parseDouble(cmd.getOptionValue("s"));
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
