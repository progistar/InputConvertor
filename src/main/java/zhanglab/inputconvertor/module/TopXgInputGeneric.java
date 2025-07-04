package zhanglab.inputconvertor.module;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public abstract class TopXgInputGeneric implements TopXgInput {
	public String inputFilePath = null;
	public String outputFilePath = null;
	public String spectrumFilePath = null;
	public String ptmTableFilePath = null;
	
	public Hashtable<String, String> loadPTMTable () throws IOException {
		Hashtable<String, String> ptmTable = new Hashtable<String, String>();
		File file = new File(ptmTableFilePath);
		
		BufferedReader BR = new BufferedReader(new FileReader(file));;
		String line = null;
		
		BR.readLine();
		// fields: 
		// Search   Unimod   ProForma
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String search = fields[0];
			String unimod = "["+fields[1]+"]";
			ptmTable.put(search, unimod);
		}
		
		BR.close();
		
		return ptmTable;
	}
	

	public void parseOptions (String[] args) {
		CommandLine cmd = null;
		Options options = new Options();
		
		// Mandatory
		Option optionInput = Option.builder("i")
				.longOpt("input").argName("res")
				.hasArg()
				.required(true)
				.desc("casanovo output file")
				.build();
		
		Option optionSpectrumFile = Option.builder("s")
				.longOpt("spectrum").argName("mgf/mzml")
				.hasArg()
				.required(true)
				.desc("spectrum file")
				.build();
		
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("file path")
				.hasArg()
				.required(true)
				.desc("output file path")
				.build();
		
		Option optionPTMTable = Option.builder("p")
				.longOpt("ptm_table").argName("tsv")
				.hasArg()
				.required(true)
				.desc("tsv file including PTM naming convention")
				.build();
		
		options.addOption(optionInput)
		.addOption(optionSpectrumFile)
		.addOption(optionOutput)
		.addOption(optionPTMTable);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
	    ArrayList<String> tmpArgs = new ArrayList<String>();
	    for(int i=0; i<args.length; i++) {
	    	if(args[i].equalsIgnoreCase("-i") || args[i].equalsIgnoreCase("--input") ||
			args[i].equalsIgnoreCase("-s") || args[i].equalsIgnoreCase("--spectrum") ||
			args[i].equalsIgnoreCase("-o") || args[i].equalsIgnoreCase("--output") ||
			args[i].equalsIgnoreCase("-p") || args[i].equalsIgnoreCase("--ptm_table")) {
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
		    
		    if(cmd.hasOption("s")) {
		    	spectrumFilePath = cmd.getOptionValue("s");
		    }
		    if(cmd.hasOption("p")) {
		    	ptmTableFilePath = cmd.getOptionValue("p");
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
