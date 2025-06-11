package zhanglab.inputconvertor.run;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.data.FastaEntry;
import zhanglab.inputconvertor.data.FastaLoader;
import zhanglab.inputconvertor.env.InputConvertorConstants;

public class RunRefineFasta {
	public static File[] referenceFiles = null;
	public static File[] filterFiles = null;
	public static File[] nonreferenceFiles = null;
	public static File outputFile = null;
	public static final Pattern PE_PATTERN = Pattern.compile("(PE=[0-9]+)");
	
	public static void main(String[] args) throws IOException, ParseException {
		System.out.println(InputConvertorConstants.VERSION);
		parseOptions(args);
		
		long startTime = System.currentTimeMillis();
		
		FastaLoader[] referenceFastas = new FastaLoader[referenceFiles.length];
		FastaLoader[] nonreferenceFastas = new FastaLoader[nonreferenceFiles.length];
		FastaLoader[] filterFastas = null;
		
		if(filterFiles != null) {
			filterFastas = new FastaLoader[filterFiles.length];
			// load filter files
			for(int i=0; i<filterFastas.length; i++) {
				filterFastas[i] = new FastaLoader(filterFiles[i]);
			}
		}
		
		System.out.println("Filter overlapped peptide sequences between the non-reference and the filter list");
		// load non-reference files
		for(int i=0; i<nonreferenceFastas.length; i++) {
			nonreferenceFastas[i] = new FastaLoader(nonreferenceFiles[i]);
			System.out.println(nonreferenceFiles[i].getName());
			// remove reference
			if(filterFiles != null) {
				for(int j=0; j<filterFastas.length; j++) {
					nonreferenceFastas[i].removeSequenceOverlapped(filterFastas[j].entries, 20);
				}
			}
		}
		
		// write database
		System.out.println("Write a final database");
		BufferedWriter BW = new BufferedWriter(new FileWriter(outputFile));
		for(int i=0; i<referenceFastas.length; i++) {
			referenceFastas[i] = new FastaLoader(referenceFiles[i]);
			
			ArrayList<FastaEntry> entries = referenceFastas[i].entries;
			for(FastaEntry entry : entries) {
				String header = entry.originHeader;
				String sequence = entry.sequence;
				
				Matcher matcher = PE_PATTERN.matcher(header);
				if(matcher.find()) {
					header = header.replace(matcher.group(), "PE=1");
				} else {
					header += " PE=1";
				}
				
				BW.append(">").append(header);
				BW.newLine();
				BW.append(sequence);
				BW.newLine();
			}
		}
		
		for(int i=0; i<nonreferenceFastas.length; i++) {
			ArrayList<FastaEntry> entries = nonreferenceFastas[i].entries;
			for(FastaEntry entry : entries) {
				String header = entry.originHeader;
				String sequence = entry.sequence;
				
				Matcher matcher = PE_PATTERN.matcher(header);
				if(matcher.find()) {
					header = header.replace(matcher.group(), "PE=2");
				} else {
					header += " PE=2";
				}
				
				BW.append(">").append(header);
				BW.newLine();
				BW.append(sequence);
				BW.newLine();
			}
			
		}
		BW.close();
        long endTime = System.currentTimeMillis();
        
        System.out.println((endTime-startTime)/1000+" sec");
	}
	
	
	public static void parseOptions (String[] args) {
		CommandLine cmd = null;
		Options options = new Options();
		
		// Mandatory
		Option optionReference = Option.builder("r")
				.longOpt("reference").argName("fasta")
				.hasArg()
				.required(true)
				.desc("a lits of reference fasta files (separted by comma). These records will be annotated as PE=1.")
				.build();
		
		Option optionFilter = Option.builder("f")
				.longOpt("filter").argName("fasta")
				.hasArg()
				.required(false)
				.desc("a lits of protein sequences to be filtered in non-reference fasta files (separted by comma).")
				.build();
		
		Option optionNonreference = Option.builder("n")
				.longOpt("non-reference").argName("fasta")
				.hasArg()
				.required(true)
				.desc("a lits of protein sequences to be filtered in non-reference fasta files (separted by comma). These records will be annotated as PE=2.")
				.build();
		
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("fasta")
				.hasArg()
				.required(true)
				.desc("a list of protein sequences after filter and merge.")
				.build();
		
		
		
		options.addOption(optionReference)
		.addOption(optionNonreference)
		.addOption(optionOutput)
		.addOption(optionFilter);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
	    ArrayList<String> tmpArgs = new ArrayList<String>();
	    for(int i=0; i<args.length; i++) {
	    	if(args[i].equalsIgnoreCase("-r") || args[i].equalsIgnoreCase("--reference") ||
			args[i].equalsIgnoreCase("-f") || args[i].equalsIgnoreCase("--filter") ||
			args[i].equalsIgnoreCase("-n") || args[i].equalsIgnoreCase("--non-reference") ||
			args[i].equalsIgnoreCase("-o") || args[i].equalsIgnoreCase("--output")) {
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
		    
		    if(cmd.hasOption("r")) {
		    	System.out.println("Reference fasta");
		    	String tmp = cmd.getOptionValue("r");
		    	String[] fileNames = tmp.split("\\,");
		    	referenceFiles = new File[fileNames.length];
		    	for(int i=0; i<referenceFiles.length; i++) {
		    		referenceFiles[i] = new File(fileNames[i]);
		    		System.out.println("  "+referenceFiles[i].getName());
		    	}
		    }
		    
		    
		    if(cmd.hasOption("n")) {
		    	System.out.println("Non-reference fasta");
		    	String tmp = cmd.getOptionValue("n");
		    	String[] fileNames = tmp.split("\\,");
		    	nonreferenceFiles = new File[fileNames.length];
		    	for(int i=0; i<nonreferenceFiles.length; i++) {
		    		nonreferenceFiles[i] = new File(fileNames[i]);
		    		System.out.println("  "+nonreferenceFiles[i].getName());
		    	}
		    }
		    

		    if(cmd.hasOption("f")) {
		    	System.out.println("Filters");
		    	String tmp = cmd.getOptionValue("f");
		    	String[] fileNames = tmp.split("\\,");
		    	filterFiles = new File[fileNames.length];
		    	for(int i=0; i<filterFiles.length; i++) {
		    		filterFiles[i] = new File(fileNames[i]);
		    		System.out.println("  "+filterFiles[i].getName());
		    	}
		    }
		    
		    if(cmd.hasOption("o")) {
		    	outputFile = new File(cmd.getOptionValue("o"));
		    	System.out.println("Output file: "+outputFile.getName());
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
