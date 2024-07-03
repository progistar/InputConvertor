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

import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.data.Annotation;
import zhanglab.inputconvertor.env.InputConvertorConstants;

public class AnnotatePeptide {

	public static ArrayList<String> inputPaths = new ArrayList<String>();
	public static String outputFilePath;
	public static String contaminantFilePath;
	public static String fdrTSVFilePath;
	
	private AnnotatePeptide () {}
	
	public static void merge(String[] args) throws IOException {
		parseOptions(args);
		
		System.out.println("Annotate peptides");
		System.out.println("- list of databases");
		for(int i=0; i<inputPaths.size(); i++) {
			System.out.println(inputPaths.get(i));
		}
		System.out.println(contaminantFilePath);
		
		File tsvFile = new File(fdrTSVFilePath);
		BufferedReader BR = new BufferedReader(new FileReader(tsvFile));
		String line = null;
		
		String header = BR.readLine(); 
		int infPeptIdx = InputConvertorConstants.getFieldIndex(header.split("\t"), InputConvertorConstants.IC_INFERRED_PEPTIDE_FIELD_NAME);
		
		ArrayList<String> records = new ArrayList<String>();
		ArrayList<String> peptides = new ArrayList<String>();
		Hashtable<String, Boolean> checkDup = new Hashtable<String, Boolean>();
		
		while((line = BR.readLine()) != null) {
			records.add(line);
			String[] fields = line.split("\t");
			String peptide = fields[infPeptIdx].replace("I", "L");
			
			if(checkDup.get(peptide) == null) {
				checkDup.put(peptide, true);
				peptides.add(peptide);
			}
		}
		
		BR.close();
		
		// information for annotated peptides
		Annotation annotation = new Annotation();
		Trie trie = Trie.builder().addKeywords(peptides).build();
		
		
		for(String fileName : inputPaths) {
			File file = new File(fileName);
			if(file.getName().startsWith(".")) continue;
			match(file, trie, annotation, false);
		}
		
		// contaminant
		File contaminantFile = new File(contaminantFilePath);
		match(contaminantFile, trie, annotation, true);

		// cal representation
		annotation.calRepresentAnnotation();
		
		File outputFile = new File(outputFilePath);
		BufferedWriter BW = new BufferedWriter(new FileWriter(outputFile));
		
		header += "\tGenomicLoci\tGenomicLociCount\tGeneName\tGeneNameCount\tMutationCount\tCategory";
		BW.append(header);
		BW.newLine();
		for(int i=0; i<records.size(); i++) {
			String record = records.get(i);
			String[] fields = record.split("\t");
			ArrayList<String> modifiedRecordList = annotation.getModifiedRecordList(record, fields[infPeptIdx]);
			
			for(String modifiedRecord : modifiedRecordList) {
				BW.append(modifiedRecord);
				BW.newLine();
			}
		}
		
		BW.close();
	}
	
	public static void match(File file, Trie trie, Annotation annotation, boolean isContam) throws IOException {

		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		String header = null;
		StringBuilder str = new StringBuilder();
		while((line = BR.readLine()) != null) {
			if(line.startsWith(">")) {
				if(header != null) {
					String sequence = str.toString();
					Collection<Emit> emits = trie.parseText(sequence.replace("I", "L"));
					for(Emit emit : emits) {
						annotation.putAnnotation(emit, header, sequence, isContam);
					}
					str.setLength(0);
				}
				
				header = line;
			} else {
				str.append(line);
			}
		}
		BR.close();
		
		// last pang
		String sequence = str.toString();
		Collection<Emit> emits = trie.parseText(sequence.replace("I", "L"));
		for(Emit emit : emits) {
			annotation.putAnnotation(emit, header, sequence, isContam);
		}
		str.setLength(0);
	}

	
	public static void parseOptions (String[] args) {
		CommandLine cmd = null;
		Options options = new Options();
		
		// Mandatory
		Option optionInput = Option.builder("i")
				.longOpt("input").argName("string")
				.hasArg()
				.required(true)
				.desc("TSV file after FDR estimation")
				.build();
		
		Option optionDatabase = Option.builder("f")
				.longOpt("fasta").argName("fa|fasta")
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
		.addOption(optionDatabase)
		.addOption(optionContaminant);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
	    ArrayList<String> tmpArgs = new ArrayList<String>();
	    for(int i=0; i<args.length; i++) {
	    	if(args[i].equalsIgnoreCase("-i") || args[i].equalsIgnoreCase("--input") ||
			args[i].equalsIgnoreCase("-o") || args[i].equalsIgnoreCase("--output") ||
			args[i].equalsIgnoreCase("-c") || args[i].equalsIgnoreCase("--contaminant") ||
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
		    cmd = parser.parse(options, nArgs, true);
		    
		    
		    if(cmd.hasOption("o")) {
		    	outputFilePath = cmd.getOptionValue("o");
		    }
		    
		    if(cmd.hasOption("i")) {
		    	fdrTSVFilePath = cmd.getOptionValue("i");
		    }
		    
		    if(cmd.hasOption("f")) {
		    	String[] inputs = cmd.getOptionValue("f").split("\\,");
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
