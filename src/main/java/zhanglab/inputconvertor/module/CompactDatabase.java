package zhanglab.inputconvertor.module;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.env.InputConvertorConstants;

public class CompactDatabase {

	public static ArrayList<String> inputPaths = new ArrayList<String>();
	public static String outputFilePath;
	public static String contaminantFilePath;
	public static String decoyPrefix = "XXX";

	private CompactDatabase () {}
	
	public static void merge(String[] args) throws IOException {
		parseOptions(args);
		
		System.out.println("Compact database");
		System.out.println("- list of databases");
		for(int i=0; i<inputPaths.size(); i++) {
			System.out.println(inputPaths.get(i));
		}
		System.out.println(contaminantFilePath);
		System.out.println("decoy prefix: "+decoyPrefix);
		
		Hashtable<String, Boolean> checkDup = new Hashtable<String, Boolean>();
		
		File outputFile = new File(outputFilePath);
		
		BufferedWriter BW = new BufferedWriter(new FileWriter(outputFile));
		
		long totalSequence = 0;
		for(String fileName : inputPaths) {
			File file = new File(fileName);
			if(file.getName().startsWith(".")) continue;
			
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String line = null;
			
			String mark = "";
			while((line = BR.readLine()) != null) {
				if(line.startsWith(">")) {
					mark = line.substring(1).split("\\_")[0];
				} else {
					totalSequence++;
					String sequence = line;
					boolean isReference = mark.startsWith(InputConvertorConstants.REF_HEADER_ID) ? true : false;
					if(checkDup.get(sequence) == null || checkDup.get(sequence) == false) {
						checkDup.put(sequence, isReference);
					}
				}
			}
			BR.close();
		}
		
		// Contaminant DB
		// Mark contaminant sequence as a reference so that it can be classified to a reference peptide in TD competition.
		BufferedReader BR = new BufferedReader(new FileReader(contaminantFilePath));
		String line = null;
		StringBuilder str = new StringBuilder();
		while((line = BR.readLine()) != null) {
			if(line.startsWith(">")) {
				if(str.length() != 0) {
					String sequence = str.toString();
					if(checkDup.get(sequence) == null) {
						checkDup.put(sequence, true);
					}
					str.setLength(0);
				}
			} else {
				str.append(line);
			}
		}
		
		BR.close();
		
		Iterator<String> sequences = (Iterator<String>) checkDup.keys();
		// target
		int idx = 1;
		while(sequences.hasNext()) {
			String sequence = sequences.next();
			String mark = InputConvertorConstants.REF_HEADER_ID;
			
			if(checkDup.get(sequence) == false) {
				mark = InputConvertorConstants.NON_REF_HEADER_ID;
			}
			
			BW.append(">"+mark+"_"+idx++);
			BW.newLine();
			BW.append(sequence);
			BW.newLine();
		}
		
		// decoy
		sequences = (Iterator<String>) checkDup.keys();
		idx = 1;
		StringBuilder sb = new StringBuilder();
		while(sequences.hasNext()) {
			String sequence = sequences.next();
			String mark = InputConvertorConstants.REF_HEADER_ID;
			
			if(checkDup.get(sequence) == false) {
				mark = InputConvertorConstants.NON_REF_HEADER_ID;
			}
			
			BW.append(">"+decoyPrefix+"_"+mark+"_"+idx++);
			BW.newLine();
			BW.append(sb.append(sequence).reverse().toString());
			BW.newLine();
			sb.setLength(0);
		}
		
		System.out.println("Compaction DB: "+totalSequence +" >> "+checkDup.size());
		
		BW.close();
	}

	
	public static void parseOptions (String[] args) {
		CommandLine cmd = null;
		Options options = new Options();
		
		// Mandatory
		Option optionInput = Option.builder("i")
				.longOpt("input").argName("string")
				.hasArg()
				.required(true)
				.desc("input fasta files to be merged (each file must be separated by comma)")
				.build();
		
		Option optionContaminant = Option.builder("c")
				.longOpt("contaminant").argName("string")
				.hasArg()
				.required(true)
				.desc("input fasta file for a contaminant database")
				.build();
		
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("fasta")
				.hasArg()
				.required(true)
				.desc("compact output fasta file")
				.build();
		
		
		options.addOption(optionInput)
		.addOption(optionOutput)
		.addOption(optionContaminant);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
	    ArrayList<String> tmpArgs = new ArrayList<String>();
	    for(int i=0; i<args.length; i++) {
	    	if(args[i].equalsIgnoreCase("-i") || args[i].equalsIgnoreCase("--input") ||
			args[i].equalsIgnoreCase("-o") || args[i].equalsIgnoreCase("--output") ||
			args[i].equalsIgnoreCase("-c") || args[i].equalsIgnoreCase("--contaminant")) {
	    		tmpArgs.add(args[i++]);
	    		tmpArgs.add(args[i]);
	    	}
	    }
	    
	    String[] nArgs = new String[tmpArgs.size()];
	    for(int i =0; i<tmpArgs.size(); i++) {
	    	nArgs[i] = tmpArgs.get(i);
	    }
	    
	    
		try {
		    cmd = parser.parse(options, nArgs, true);
		    
		    
		    if(cmd.hasOption("o")) {
		    	outputFilePath = cmd.getOptionValue("o");
		    }
		    
		    if(cmd.hasOption("i")) {
		    	String[] inputs = cmd.getOptionValue("i").split("\\,");
		    	for(String input : inputs) {
		    		inputPaths.add(input);
		    	}
		    }
		    
		    if(cmd.hasOption("c")) {
		    	contaminantFilePath = cmd.getOptionValue("c");
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
