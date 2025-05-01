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

public class ToModModifier {

	public static String inputFilePath;
	public static String outputFilePath;
	public static String ptmTableFilePath;
	
	private ToModModifier() {};
	
	public static void modifyPTM (String[] args) throws IOException {
		parseOptions(args);
		Hashtable<String, String> ptmTable = loadPTMTable();
		ArrayList<String> ptms = new ArrayList<String>();
		ptmTable.forEach((s, u)->{
			ptms.add(s);
		});
		
		BufferedReader BR = new BufferedReader(new FileReader(inputFilePath));
		BufferedWriter BW = new BufferedWriter(new FileWriter(outputFilePath));
		String line = null;
		
		String header = BR.readLine();
		
		String[] fields = header.split("\t");
		int peptideIdx = -1;
		
		for(int i=0; i<fields.length; i++) {
			if(fields[i].equalsIgnoreCase("Peptide")) {
				peptideIdx = i;
			}
		}
		
		BW.append(header);
		BW.newLine();
		while((line = BR.readLine()) != null) {
			fields = line.split("\t");
			
			// replace all PTMs
			for(String ptm : ptms) {
				fields[peptideIdx] = fields[peptideIdx].replace(ptm, ptmTable.get(ptm));
			}
			
			BW.append(fields[0]);
			for(int i=1; i<fields.length; i++) {
				BW.append("\t").append(fields[i]);
			}
			BW.newLine();
		}
		
		
		BW.close();
		BR.close();
	}
	
	public static Hashtable<String, String> loadPTMTable () throws IOException {
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
	
	public static void parseOptions (String[] args) {
		CommandLine cmd = null;
		Options options = new Options();
		
		// Mandatory
		Option optionInputTxt = Option.builder("i")
				.longOpt("input").argName("pin")
				.hasArg()
				.required(true)
				.desc("feat file")
				.build();
		
		Option optionPTMTable = Option.builder("p")
				.longOpt("ptm_table").argName("tsv")
				.hasArg()
				.required(true)
				.desc("tsv file including PTM naming convention")
				.build();
		
		
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("pin")
				.hasArg()
				.required(true)
				.desc("pin output file")
				.build();
		options.addOption(optionInputTxt)
		.addOption(optionPTMTable)
		.addOption(optionOutput);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
	    ArrayList<String> tmpArgs = new ArrayList<String>();
	    for(int i=0; i<args.length; i++) {
	    	if(args[i].equalsIgnoreCase("-i") || args[i].equalsIgnoreCase("--input") ||
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
