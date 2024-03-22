package zhanglab.inputconvertor.module;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
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

import zhanglab.inputconvertor.env.InputConvertorConstants;

public class ToNetMHCpanInput {

	public static String inputFilePath;
	public static String outputFilePath;
	
	private ToNetMHCpanInput () {};
	public static void makeNetMHCpanInput (String[] args) throws IOException {
		parseOptions(args);
		
		int infPeptideIdx		= -1;
		
		File inputFile = new File(inputFilePath);
		
		System.out.println("Make NetMHCpanInput file from "+inputFile.getName());
		
		BufferedReader BR = new BufferedReader(new FileReader(inputFile));
		BufferedWriter BW = new BufferedWriter(new FileWriter(outputFilePath));
		
		String line = null;
		String[] header = BR.readLine().split("\t");
		
		infPeptideIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.IC_INFERRED_PEPTIDE_FIELD_NAME);
		
		Hashtable<String, String> isDuplicated = new Hashtable<String, String>();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String infPeptide = fields[infPeptideIdx];
			
			if(infPeptide.length() >= 8 && infPeptide.length() <= 15) {

    			if(isDuplicated.get(infPeptide) == null) {
    				BW.append(infPeptide);
    				BW.newLine();
    				isDuplicated.put(infPeptide, "");
    			}
			}
			
		}
		
		System.out.println("A total of "+isDuplicated.size()+" peptides were collected from the file.");
		
		BR.close();
		BW.close();
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
				.desc("output peptide list file")
				.build();
		
		
		options.addOption(optionInput)
		.addOption(optionOutput);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
	    ArrayList<String> tmpArgs = new ArrayList<String>();
	    for(int i=0; i<args.length; i++) {
	    	if(args[i].equalsIgnoreCase("-i") || args[i].equalsIgnoreCase("--input") ||
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
		    cmd = parser.parse(options, nArgs, false);
		    
		    if(cmd.hasOption("i")) {
		    	inputFilePath = cmd.getOptionValue("i");
		    }
		    
		    if(cmd.hasOption("o")) {
		    	outputFilePath = cmd.getOptionValue("o");
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

