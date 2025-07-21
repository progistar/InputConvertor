package zhanglab.inputconvertor.run;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Pattern;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.data.FragPipeWorkflow;
import zhanglab.inputconvertor.env.InputConvertorConstants;

public class RunMakeFragpipeInput {
	
	public static String sampleName = null;
	public static File workflowFile = null;
	public static File[] spectrumFiles = null;
	public static File neoDBFile = null;
	public static String msType = null;
	public static File outputWorkflowFile = null;
	public static File outputManifestFile = null;
	public static String decoyPrefix = null;
	public static final Pattern PE_PATTERN = Pattern.compile("(PE=[0-9]+)");
	
	public static void main(String[] args) throws IOException, ParseException {
		System.out.println(InputConvertorConstants.VERSION);
		parseOptions(args);
		
		long startTime = System.currentTimeMillis();
		
		// read workflow file
		FragPipeWorkflow workflow = new FragPipeWorkflow(workflowFile);
		int sliceNum = estimateSliceDBNum(neoDBFile);
		
		// set params
		workflow.setParam("database.db-path", "/data/"+neoDBFile.getName()); // set db path
		workflow.setParam("msfragger.misc.slice-db", sliceNum+"");  // set slice-db
		workflow.setParam("database.decoy-tag", decoyPrefix); // set decoy prefix
		workflow.setParam("phi-report.print-decoys", "true"); // print decoy
		workflow.setParam("msfragger.group_variable", "2"); // set separate FDR
		
		
		// write new workflow file
		workflow.write(outputWorkflowFile);
		
		// write new manifest file
		writeManifest(outputManifestFile, spectrumFiles, msType, sampleName);
		
        long endTime = System.currentTimeMillis();
        
        System.out.println((endTime-startTime)/1000+" sec");
	}
	
	public static void writeManifest (
			File outputManifestFile,
			File[] spectrumFiles,
			String msType,
			String sampleName) throws IOException {
		
		BufferedWriter BW = new BufferedWriter(new FileWriter(outputManifestFile));
		
		for(File file : spectrumFiles) {
			BW.append("/data/"+file.getName()).append("\t").append(sampleName).append("\t\t").append(msType);
			BW.newLine();
		}
		
		BW.close();
		
	}
	
	public static int estimateSliceDBNum (File neoDBFile) throws IOException {
		int slice = 1;
		long aaCnt = 0;
		
		BufferedReader BR = new BufferedReader(new FileReader(neoDBFile));
		String line = null;
		
		while((line = BR.readLine()) != null) {
			if(line.startsWith(">")) continue;
			aaCnt += line.length();
		}
		
		BR.close();
		
		slice = (int)(aaCnt/100000000) + 1;
		
		System.out.println("Total AA counts: "+aaCnt);
		System.out.println("Estimated slice num: "+slice);
		
		return slice;
	}
	
	public static void parseOptions (String[] args) {
		CommandLine cmd = null;
		Options options = new Options();
		
		// Mandatory
		Option optionWorkflow = Option.builder("w")
				.longOpt("workflow").argName("workflow")
				.hasArg()
				.required(true)
				.desc("FragPipe workflow file.")
				.build();
		
		Option optionDecoy = Option.builder("d")
				.longOpt("decoy").argName("string")
				.hasArg()
				.required(true)
				.desc("Specify decoy prefix. Decoy sequences are generated only if the value is specified. For example, you can set rev_ for FragPipe.")
				.build();
		
		Option optionSpectrumPath = Option.builder("s")
				.longOpt("spectrum").argName("string")
				.hasArg()
				.required(true)
				.desc("A list of spectrum files separated by comma.")
				.build();
		
		Option optionNeoDBPath = Option.builder("db")
				.longOpt("neodb").argName("string")
				.hasArg()
				.required(true)
				.desc("Path of neodb.")
				.build();
		
		Option optionSampleName = Option.builder("n")
				.longOpt("name").argName("string")
				.hasArg()
				.required(true)
				.desc("Sample name.")
				.build();
		
		Option optionMSType = Option.builder("m")
				.longOpt("ms").argName("string")
				.hasArg()
				.required(true)
				.desc("A type of MS.")
				.build();
		
		
		Option optionOutputWorkflow = Option.builder("ow")
				.longOpt("output_workflow").argName("string")
				.hasArg()
				.required(true)
				.desc("Output path of new workflow.")
				.build();
		
		Option optionOutputManifest = Option.builder("om")
				.longOpt("output_manifest").argName("string")
				.hasArg()
				.required(true)
				.desc("Output path of new manifest.")
				.build();
		
		
		options.addOption(optionWorkflow)
		.addOption(optionSpectrumPath)
		.addOption(optionSampleName)
		.addOption(optionDecoy)
		.addOption(optionMSType)
		.addOption(optionNeoDBPath)
		.addOption(optionOutputWorkflow)
		.addOption(optionOutputManifest);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
	    ArrayList<String> tmpArgs = new ArrayList<String>();
	    for(int i=0; i<args.length; i++) {
	    	if(args[i].equalsIgnoreCase("-w") || args[i].equalsIgnoreCase("--workflow") ||
			args[i].equalsIgnoreCase("-s") || args[i].equalsIgnoreCase("--spectrum") ||
			args[i].equalsIgnoreCase("-n") || args[i].equalsIgnoreCase("--name") ||
			args[i].equalsIgnoreCase("-d") || args[i].equalsIgnoreCase("--decoy") ||
			args[i].equalsIgnoreCase("-db") || args[i].equalsIgnoreCase("--neodb") ||
			args[i].equalsIgnoreCase("-m") || args[i].equalsIgnoreCase("--ms") ||
			args[i].equalsIgnoreCase("-ow") || args[i].equalsIgnoreCase("--output_workflow") ||
			args[i].equalsIgnoreCase("-om") || args[i].equalsIgnoreCase("--output_manifest")) {
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
		    
		    if(cmd.hasOption("n")) {
		    	System.out.println("Sample name");
		    	String tmp = cmd.getOptionValue("n");
		    	sampleName = tmp;
		    	System.out.println("  "+sampleName);
		    }
		    
		    if(cmd.hasOption("w")) {
		    	System.out.println("Workflow file");
		    	String tmp = cmd.getOptionValue("w");
		    	workflowFile = new File(tmp);
		    	System.out.println("  "+workflowFile.getName());
		    }
		    

		    if(cmd.hasOption("s")) {
		    	System.out.println("Spectrum files");
		    	String tmp = cmd.getOptionValue("s");
		    	String[] fileNames = tmp.split("\\,");
		    	spectrumFiles = new File[fileNames.length];
		    	for(int i=0; i<spectrumFiles.length; i++) {
		    		spectrumFiles[i] = new File(fileNames[i].trim());
		    		System.out.println("  "+spectrumFiles[i].getName());
		    	}
		    }
		    
		    if(cmd.hasOption("m")) {
		    	System.out.println("MS type");
		    	msType = cmd.getOptionValue("m");
		    	System.out.println("  "+msType);
		    }
		    
		    
		    if(cmd.hasOption("d")) {
		    	System.out.println("Decoy prefix");
		    	decoyPrefix = cmd.getOptionValue("d");
		    	System.out.println("  "+decoyPrefix);
		    }
		    
		    if(cmd.hasOption("db")) {
		    	System.out.println("NeoDB");
		    	String tmp = cmd.getOptionValue("db");
		    	neoDBFile = new File(tmp);
		    	System.out.println("  "+neoDBFile.getName());
		    }
		    
		    if(cmd.hasOption("ow")) {
		    	System.out.println("Output workflow file");
		    	String tmp = cmd.getOptionValue("ow");
		    	outputWorkflowFile = new File(tmp);
		    	System.out.println("  "+outputWorkflowFile.getName());
		    }
		    
		    if(cmd.hasOption("om")) {
		    	System.out.println("Output manifest file");
		    	String tmp = cmd.getOptionValue("om");
		    	outputManifestFile = new File(tmp);
		    	System.out.println("  "+outputManifestFile.getName());
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
