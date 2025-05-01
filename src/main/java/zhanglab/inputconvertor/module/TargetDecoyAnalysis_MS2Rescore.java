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

public class TargetDecoyAnalysis_MS2Rescore {

	public static String inputFilePath = null;
	public static String outputFilePath = null;
	public static String ms2RescoreFilePath = null;
	public static String pinFilePath = null;
	public static double fdr = -1;
	public static boolean isLengthSpecific = false;
	public static boolean isPrintDecoy = false;
	
	private TargetDecoyAnalysis_MS2Rescore() {};
	
	
	public static void doFDR (String[] args) throws IOException, ParseException {
		
        parseOptions(args);
        
        Hashtable<String, String[]> topRankPSMs = readMS2RescoreTSV(ms2RescoreFilePath);
        
        
        BufferedReader BR = new BufferedReader(new FileReader(inputFilePath));
        BufferedReader PIN = new BufferedReader(new FileReader(pinFilePath));
        BufferedWriter BW = new BufferedWriter(new FileWriter(outputFilePath));
        String line = null;
        String header = BR.readLine();
        String[] pinHeader = PIN.readLine().split("\t");
        
        int psmIdIdx = InputConvertorConstants.getFieldIndex(pinHeader, "SpecId");
		int proteinIdIdx = InputConvertorConstants.getFieldIndex(pinHeader, "Proteins");
		int peptideIdx = InputConvertorConstants.getFieldIndex(pinHeader, "Peptide");
        
        // build header ////////////////////////////////////////////////
        BW.append(header).append("\t")
        .append("final_score");
        BW.newLine();
        ////////////////////////////////////// End of building header /
        
        String[] headerSplit = header.split("\t");
        int specIdIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.IC_SPEC_ID_FEILD_NAME);
        int labelIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.IC_LABEL_FEILD_NAME);
        int isCanonicalIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.IC_IS_CANONICAL_FEILD_NAME);
        int inferredPeptideIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.IC_INFERRED_PEPTIDE_FIELD_NAME);

        
        ArrayList<TDRecord> cRecords = new ArrayList<TDRecord>();
        ArrayList<TDRecord> ncRecords = new ArrayList<TDRecord>();
        
        Hashtable<String, String> duplications = new Hashtable<String, String>();
        while((line = BR.readLine()) != null) {
        	String pinRecord = PIN.readLine();
        	
        	// generate key from pin file
        	String[] fields = pinRecord.split("\t");
    		
    		StringBuilder proteins = new StringBuilder("[");
    		for(int i=proteinIdIdx; i<fields.length; i++) {
    			proteins.append("'").append(fields[i]).append("',");
    		}
    		proteins.setLength(proteins.length()-1); // remove , at the end of string
    		proteins.append("]");
    		
    		String key = fields[psmIdIdx]+"+"+fields[peptideIdx]+"+"+proteins.toString();
        	
    		String[] ms2RescoreRecord = topRankPSMs.get(key);
    		if(ms2RescoreRecord == null) {
    			continue;
    		}
    		
    		
    		
        	fields = line.split("\t");
        	String isCanonical = fields[isCanonicalIdx];
        	String label = fields[labelIdx];
        	int peptideLength = fields[inferredPeptideIdx].replace("[+-0123456789.*\\(\\)]", "").length();
        	TDRecord record = new TDRecord();
    		record.record = line;
    		record.score = Double.parseDouble(ms2RescoreRecord[0]);
    		record.label = Integer.parseInt(label);
    		record.peptideLength = peptideLength;
    		
    		if(isCanonical.equalsIgnoreCase("true")) {
    			cRecords.add(record);
    		} else {
    			ncRecords.add(record);
    		}
        }
        
        // length-specific FDR control
        ArrayList<TDRecord> passList = new ArrayList<TDRecord>();
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
            	passList.addAll(getFDR(lengthSpecificCRecords, fdr));
            	System.out.println("Non-canonical PSMs with length "+len);
            	System.out.println("Target PSMs: "+getNumOfRecords(lengthSpecificNCRecords, 1));
            	System.out.println("Decoy PSMs: "+getNumOfRecords(lengthSpecificNCRecords, -1));
            	passList.addAll(getFDR(lengthSpecificNCRecords, fdr));
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
        	passList.addAll(pass);
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
        	passList.addAll(pass);
        }
        
        
        for(int i=0; i<passList.size(); i++) {
        	TDRecord record = passList.get(i);
        	
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
        PIN.close();
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
		
		ArrayList<TDRecord> passList = new ArrayList<TDRecord>();
		for(int i=0; i<=lastIdx; i++) {
			TDRecord record = records.get(i);
			if(isPrintDecoy) {
				passList.add(record);
			} else if(record.label > 0) {
				passList.add(record);
			}
		}
		
		System.out.println("FDR: " + estimatedFDR);
		System.out.println("IDs: "+passList.size());
		
		return passList;
	}
	
	
	
	private static Hashtable<String, String[]> readMS2RescoreTSV (String fileName) throws IOException {
		File file = new File(fileName);
		Hashtable<String, String[]> table = new Hashtable<String, String[]>();
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		String[] header = BR.readLine().split("\t");
		int psmIdIdx = InputConvertorConstants.getFieldIndex(header, "spectrum_id");
		int scoreIdx = InputConvertorConstants.getFieldIndex(header, "score");
		int proteinIdIdx = InputConvertorConstants.getFieldIndex(header, "protein_list");
		int peptideIdx = InputConvertorConstants.getFieldIndex(header, "peptidoform");
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String[] record = {fields[scoreIdx], line};
			// peptidoform = Peptide/Charge
			String peptide = fields[peptideIdx].split("\\/")[0];
			String key = fields[psmIdIdx]+"+"+peptide+"+"+fields[proteinIdIdx];
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
		
		Option optionPinFile = Option.builder("p")
				.longOpt("pin").argName("pin")
				.hasArg()
				.required(true)
				.desc("input pin file for MS2Rescore")
				.build();
		
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("tsv")
				.hasArg()
				.required(true)
				.desc("fdr output file")
				.build();
		
		Option optionMS2RescoreTSV = Option.builder("t")
				.longOpt("ms2rescore").argName("tsv")
				.hasArg()
				.required(true)
				.desc("MS2Rescore result file")
				.build();
		
		Option optionFDR = Option.builder("f")
				.longOpt("fdr").argName("float")
				.hasArg()
				.required(true)
				.desc("FDR value")
				.build();
		
		Option optionLengthSpecificFDR = Option.builder("l")
				.longOpt("length_specific")
				.required(false)
				.desc("Calculate a length-specific FDR")
				.build();
		
		Option optionPrintDecoy = Option.builder("d")
				.longOpt("print_decoy")
				.required(false)
				.desc("Print decoys in the final result")
				.build();
		
		options.addOption(optionInput)
		.addOption(optionPinFile)
		.addOption(optionMS2RescoreTSV)
		.addOption(optionFDR)
		.addOption(optionOutput)
		.addOption(optionPrintDecoy)
		.addOption(optionLengthSpecificFDR);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
	    ArrayList<String> tmpArgs = new ArrayList<String>();
	    for(int i=0; i<args.length; i++) {
	    	if(args[i].equalsIgnoreCase("-i") || args[i].equalsIgnoreCase("--input") ||
			args[i].equalsIgnoreCase("-t") || args[i].equalsIgnoreCase("--ms2rescore") ||
			args[i].equalsIgnoreCase("-p") || args[i].equalsIgnoreCase("--pin") ||
			args[i].equalsIgnoreCase("-f") || args[i].equalsIgnoreCase("--fdr") ||
			args[i].equalsIgnoreCase("-o") || args[i].equalsIgnoreCase("--output")) {
	    		tmpArgs.add(args[i++]);
	    		tmpArgs.add(args[i]);
	    	} else if(args[i].equalsIgnoreCase("-l") || args[i].equalsIgnoreCase("--length_specific") ||
	    			args[i].equalsIgnoreCase("-d") || args[i].equalsIgnoreCase("--print_decoy")) {
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
		    
		    if(cmd.hasOption("p")) {
		    	pinFilePath = cmd.getOptionValue("p");
		    }
		    
		    if(cmd.hasOption("t")) {
		    	ms2RescoreFilePath = cmd.getOptionValue("t");
		    }
		    
		    if(cmd.hasOption("f")) {
		    	fdr = Double.parseDouble(cmd.getOptionValue("f"));
		    }
		    
		    if(cmd.hasOption("l")) {
		    	isLengthSpecific = true;
		    }
		    if(cmd.hasOption("d")) {
		    	isPrintDecoy = true;
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
