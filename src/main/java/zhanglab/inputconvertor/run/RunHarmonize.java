package zhanglab.inputconvertor.run;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.data.FastaLoader;
import zhanglab.inputconvertor.data.HarmonyFragPipe;
import zhanglab.inputconvertor.data.HarmonypXg;
import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.input.HarmonyResult;

public class RunHarmonize {
	public static File[] inputFiles = null;
	public static File fastaFile = null;
	public static File outputFile = null;
	public static String pipeline = null;
	public static String decoyPrefix = "rev_";
	
	public static void main(String[] args) throws IOException, ParseException {
		System.out.println(InputConvertorConstants.VERSION);
		parseOptions(args);
		
		long startTime = System.currentTimeMillis();
		HarmonyResult[] harmonyResults = new HarmonyResult[inputFiles.length];
		
		for(int i=0; i<inputFiles.length; i++) {
			if(pipeline.equalsIgnoreCase(InputConvertorConstants.PXG)) {
				harmonyResults[i] = new HarmonypXg(inputFiles[i]);
			} else if(pipeline.equalsIgnoreCase(InputConvertorConstants.FRAGPIPE)) {
				harmonyResults[i] = new HarmonyFragPipe(inputFiles[i]);
			}
		}
		
		// merge
		HarmonyResult mergedResult = HarmonyResult.merge(harmonyResults);
		
		FastaLoader fLoader = new FastaLoader(fastaFile);

		mergedResult.runSequentialSteps(fLoader, decoyPrefix);
		// write a final result
		mergedResult.write(outputFile);
		
        long endTime = System.currentTimeMillis();
        
        System.out.println((endTime-startTime)/1000+" sec");
	}
	
	
	public static void parseOptions (String[] args) {
		CommandLine cmd = null;
		Options options = new Options();
		
		// Mandatory
		Option optionInput = Option.builder("i")
				.longOpt("input").argName("tsv")
				.hasArg()
				.required(true)
				.desc("list of annotations from PepQueryMHC.")
				.build();
		
		Option optionFasta = Option.builder("f")
				.longOpt("fasta").argName("fasta")
				.hasArg()
				.required(true)
				.desc("path of sequence database.")
				.build();
		
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("string")
				.hasArg()
				.required(true)
				.desc("harmonized result.")
				.build();
		
		Option optionPipeline = Option.builder("p")
				.longOpt("pipeline").argName("pXg|fragpipe")
				.hasArg()
				.required(true)
				.desc("specify a search pipeline (pXg|fragpipe).")
				.build();
		
		Option optionDecoy = Option.builder("d")
				.longOpt("decoy").argName("string")
				.hasArg()
				.required(false)
				.desc("specify decoy prefix (default is rev_).")
				.build();
		
		
		
		options.addOption(optionInput)
		.addOption(optionOutput)
		.addOption(optionDecoy)
		.addOption(optionPipeline)
		.addOption(optionFasta);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
	    ArrayList<String> tmpArgs = new ArrayList<String>();
	    for(int i=0; i<args.length; i++) {
	    	if(args[i].equalsIgnoreCase("-i") || args[i].equalsIgnoreCase("--input") ||
			args[i].equalsIgnoreCase("-f") || args[i].equalsIgnoreCase("--fasta") ||
			args[i].equalsIgnoreCase("-p") || args[i].equalsIgnoreCase("--pipeline") ||
			args[i].equalsIgnoreCase("-d") || args[i].equalsIgnoreCase("--decoy") ||
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
		    
		    if(cmd.hasOption("f")) {
		    	fastaFile = new File(cmd.getOptionValue("f"));
		    	System.out.println("Fasta file: "+fastaFile.getName());
		    }
		    
		    if(cmd.hasOption("i")) {
		    	String[] filePaths = cmd.getOptionValue("i").split("\\,");
		    	inputFiles = new File[filePaths.length];
		    	for(int i=0; i<filePaths.length; i++) {
		    		inputFiles[i] = new File(filePaths[i]);
		    		System.out.println("Input file"+(i+1)+": "+inputFiles[i].getName());
		    	}
		    	
		    }
		    
		    if(cmd.hasOption("p")) {
		    	pipeline = cmd.getOptionValue("p");
		    	if(pipeline.equalsIgnoreCase(InputConvertorConstants.PXG) ||
		    			pipeline.equalsIgnoreCase(InputConvertorConstants.FRAGPIPE)) {
		    		System.out.println("Pipeline: "+pipeline);
		    	} else {
		    		System.out.println("Use pXg or fragpipe");
		    		isFail = true;
		    	}
		    }
		    
		    if(cmd.hasOption("d")) {
		    	System.out.println("Decoy prefix");
		    	decoyPrefix = cmd.getOptionValue("d");
		    	System.out.println("  "+decoyPrefix);
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
