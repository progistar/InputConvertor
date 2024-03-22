package zhanglab.inputconvertor.module;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.data.NetMHCpanResult;
import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.function.NetMHCpanParser;

public class MergeNetMHCpan {
	
	public static String inputFilePath;
	public static String netMHCpanPath;
	public static String outputFilePath;

	private MergeNetMHCpan () {}
	
	public static void merge(String[] args) throws IOException {
		parseOptions(args);
		
		int infPeptideIdx		= -1;
		File inputFile = new File(inputFilePath);
		File netMHCFile = new File(netMHCpanPath);

        // Read NetMHCpan Output ////////////////////////////////////////////////////////
        NetMHCpanResult result = NetMHCpanParser.parseNetMHCpan(netMHCFile.getAbsolutePath());
        /////////////////////////////////// End of NetMHCpan Output reader //////////////
		
		BufferedWriter BW = new BufferedWriter(new FileWriter(outputFilePath));
		BufferedReader BR = new BufferedReader(new FileReader(inputFile));
		String line = null;
		String header = BR.readLine();
		String[] headerSplit = header.split("\t");
		
		infPeptideIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.IC_INFERRED_PEPTIDE_FIELD_NAME);
		
		BW.append(header).append("\t").append(result.getHeader());
		BW.newLine();
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String infPeptide = fields[infPeptideIdx];
			
			String ba = result.getHLATyping(infPeptide);
			
			BW.append(line).append("\t").append(ba);
			BW.newLine();
		}
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
		
		Option optionNetMHCpan = Option.builder("x")
				.longOpt("netmhcpan").argName("xls")
				.hasArg()
				.required(true)
				.desc("output file of NetMHCpan 4.1")
				.build();
		
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("tsv")
				.hasArg()
				.required(true)
				.desc("output peptide list file")
				.build();
		
		
		options.addOption(optionInput)
		.addOption(optionNetMHCpan)
		.addOption(optionOutput);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
	    ArrayList<String> tmpArgs = new ArrayList<String>();
	    for(int i=0; i<args.length; i++) {
	    	if(args[i].equalsIgnoreCase("-i") || args[i].equalsIgnoreCase("--input") ||
			args[i].equalsIgnoreCase("-o") || args[i].equalsIgnoreCase("--output") ||
			args[i].equalsIgnoreCase("-x") || args[i].equalsIgnoreCase("--netmhcpan")) {
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
		    
		    if(cmd.hasOption("x")) {
		    	netMHCpanPath = cmd.getOptionValue("x");
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
