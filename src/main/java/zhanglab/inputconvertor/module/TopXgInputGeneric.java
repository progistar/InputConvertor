package zhanglab.inputconvertor.module;

import java.util.ArrayList;

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
	public boolean isAppend = false;
	

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
		
		Option optionAppend = Option.builder("a")
				.longOpt("append")
				.required(false)
				.desc("append to exist file")
				.build();
		
		
		options.addOption(optionInput)
		.addOption(optionSpectrumFile)
		.addOption(optionOutput)
		.addOption(optionAppend);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
	    ArrayList<String> tmpArgs = new ArrayList<String>();
	    for(int i=0; i<args.length; i++) {
	    	if(args[i].equalsIgnoreCase("-i") || args[i].equalsIgnoreCase("--input") ||
			args[i].equalsIgnoreCase("-s") || args[i].equalsIgnoreCase("--spectrum") ||
			args[i].equalsIgnoreCase("-o") || args[i].equalsIgnoreCase("--output")) {
	    		tmpArgs.add(args[i++]);
	    		tmpArgs.add(args[i]);
	    	} else if (args[i].equalsIgnoreCase("-a") || args[i].equalsIgnoreCase("--append")) {
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
		    if(cmd.hasOption("a")) {
		    	isAppend = true;
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
