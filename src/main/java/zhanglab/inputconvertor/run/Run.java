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
import zhanglab.inputconvertor.input.pNovo3;
import zhanglab.inputconvertor.module.CompactDatabase;
import zhanglab.inputconvertor.module.ConvertCometTopXg;
import zhanglab.inputconvertor.module.MergeNetMHCpan;
import zhanglab.inputconvertor.module.TargetDecoyAnalysis;
import zhanglab.inputconvertor.module.ToAutoRTInput;
import zhanglab.inputconvertor.module.ToFeatures;
import zhanglab.inputconvertor.module.ToFeaturesComet;
import zhanglab.inputconvertor.module.ToNetMHCpanInput;
import zhanglab.inputconvertor.module.ToPIN;
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
        	}
        	
        	topXgInput.topXgInputFormat(args);
        }
        else if(module.equalsIgnoreCase("autort_input")) {
        	ToAutoRTInput.toAutoRTInputFormat(args, tool);
        } 
        else if(module.equalsIgnoreCase("feature_gen_denovo")) {
        	ToFeatures.toFeatureFormat(args);
        }
        else if(module.equalsIgnoreCase("feature_gen_db")) {
        	ToFeaturesComet.toFeatureFormat(args);
        }
        else if(module.equalsIgnoreCase("feature_to_pin")) {
        	ToPIN.convertToPIN(args);
        }
        else if(module.equalsIgnoreCase("fdr")) {
        	// -i -p -d
        	TargetDecoyAnalysis.doFDR(args);
        }
        else if(module.equalsIgnoreCase("netmhcpan_input")) {
        	ToNetMHCpanInput.makeNetMHCpanInput(args);
        }
        else if(module.equalsIgnoreCase("netmhcpan_merge")) {
        	MergeNetMHCpan.merge(args);
        }
        else if(module.equalsIgnoreCase("convert_to_pxg")) {
        	ConvertCometTopXg.convertTopXgOutput(args);
        }
        else if(module.equalsIgnoreCase("compact_database")) {
        	CompactDatabase.merge(args);
        }
        /*
        else if(mode.equalsIgnoreCase("ms2pip_input")) {
        	// args: -i -p -C
        	// -i ./
        	// -p C3N-01488.T.csnv
        	// --charge 14
        	ToMS2PIPInput toMS2PIPInput = null;
        	if(inputTool.equalsIgnoreCase(InputConvertorConstants.PXG)) {
        		toMS2PIPInput = new pXg();
        	} else if(inputTool.equalsIgnoreCase(InputConvertorConstants.COMET)) {
        		toMS2PIPInput = new Comet();
        	}
        	
        	toMS2PIPInput.toMS2PIPInputFormat(cmd);
        	
        }
        
        
        else if(mode.equalsIgnoreCase("fdr")) {
        	// -i -p -d
        	new TargetDecoyAnalysis().doFDR(cmd);
        }
        
        else if(mode.equalsIgnoreCase("netmhcpan_input")) {
        	NetMHCpan toNetMHCpanInput = null;
        	
        	toNetMHCpanInput = new pXg();
        	
        	toNetMHCpanInput.toNetMHCpanInputFormat(cmd);
        }
        
        else if(mode.equalsIgnoreCase("add_netmhcpan")) {
        	NetMHCpan toNetMHCpanInput = null;
        	
        	toNetMHCpanInput = new pXg();
        	
        	toNetMHCpanInput.addNetMHCpanOutput(cmd);
        }
        // this is hidden function for me
        // Do not use general purpose!!!
        else if(mode.equalsIgnoreCase("spectrum_count")) {
        	// args: -f
        	new SpectrumCount(cmd);
        }
        
        else if(mode.equalsIgnoreCase("translate")) {
        	ToTranslator toTranslator = null;
        	
        	toTranslator.doTranslate(cmd);
        }
        
        else if(mode.equalsIgnoreCase("gen_comet_param")) {
        	
        	new GenerateCometParams(args);
        }
        */
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
