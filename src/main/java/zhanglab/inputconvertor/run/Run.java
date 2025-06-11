package zhanglab.inputconvertor.run;

import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.input.Casanovo;
import zhanglab.inputconvertor.input.PEAKS;
import zhanglab.inputconvertor.input.PEAKS12;
import zhanglab.inputconvertor.input.PiPrimeNovo;
import zhanglab.inputconvertor.input.pNovo3;
import zhanglab.inputconvertor.module.TargetDecoyAnalysis_MS2Rescore;
import zhanglab.inputconvertor.module.TopXgInput;

public class Run {

	public static String module = null;
	public static String tool = null;
	
	public static void main(String[] args) throws IOException, ParseException {
		long startTime = System.currentTimeMillis();
		System.out.println(InputConvertorConstants.VERSION);
		parseOptions(args);
		
        
        if(module.equalsIgnoreCase("pxg_input")) {
    		// args: -i -p -b -f
        	TopXgInput topXgInput = null;
        	if(tool.equalsIgnoreCase(InputConvertorConstants.CASANOVO)) {
        		topXgInput = new Casanovo();
        	} else if(tool.equalsIgnoreCase(InputConvertorConstants.PNOVO3)) {
        		topXgInput = new pNovo3();
        	} else if(tool.equalsIgnoreCase(InputConvertorConstants.PEAKS)) {
        		topXgInput = new PEAKS();
        	} else if(tool.equalsIgnoreCase(InputConvertorConstants.PEAKS12)) {
        		topXgInput = new PEAKS12();
        	} else if(tool.equalsIgnoreCase(InputConvertorConstants.piPRIME_NOVO)) {
        		topXgInput = new PiPrimeNovo();
        	}
        	
        	topXgInput.topXgInputFormat(args);
        }
        else if(module.equalsIgnoreCase("ms2rescore")) {
        	// -i -p -d
        	TargetDecoyAnalysis_MS2Rescore.doFDR(args, InputConvertorConstants.PSM_LEVEL);
        	TargetDecoyAnalysis_MS2Rescore.doFDR(args, InputConvertorConstants.PEPTIDE_LEVEL);
        }
        
        long endTime = System.currentTimeMillis();
        
        System.out.println((endTime-startTime)/1000+" sec");
	}
	
	public static void parseOptions (String[] args) {
		CommandLine cmd = null;
		Options options = new Options();
		
		// Mandatory
		Option optionModule = Option.builder("m")
				.longOpt("module").argName("pxg_input/")
				.hasArg()
				.required(true)
				.desc("specify a module to use")
				.build();
		
		Option optionTool = Option.builder("t")
				.longOpt("tool").argName("casanovo/pnovo3/peaks/comet/pxg")
				.hasArg()
				.required(false)
				.desc("specify a format of result to convert")
				.build();
		
		options.addOption(optionModule)
		.addOption(optionTool);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
	    ArrayList<String> argList = new ArrayList<String>();
	    for(int i=0; i<args.length; i++) {
	    	if(args[i].equalsIgnoreCase("-m") || args[i].equalsIgnoreCase("--module") ||
			args[i].equalsIgnoreCase("-t") || args[i].equalsIgnoreCase("--tool")) {
	    		argList.add(args[i++]);
	    		argList.add(args[i]);
	    	}
	    }
	    
	    
	    String[] nArgs = new String[argList.size()];
	    for(int i=0; i<nArgs.length; i++) {
	    	nArgs[i] = argList.get(i);
	    }
	    
		try {
		    cmd = parser.parse(options, nArgs, false);
		    
		    if(cmd.hasOption("m")) {
		    	module = cmd.getOptionValue("m");
		    }
		    
		    if(cmd.hasOption("t")) {
		    	tool = cmd.getOptionValue("t");
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
