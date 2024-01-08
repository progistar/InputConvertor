package zhanglab.inputconvertor.run;

import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.input.Arriba;
import zhanglab.inputconvertor.input.Casanovo;
import zhanglab.inputconvertor.input.Comet;
import zhanglab.inputconvertor.input.PEAKS;
import zhanglab.inputconvertor.input.pNovo3;
import zhanglab.inputconvertor.input.pXg;
import zhanglab.inputconvertor.module.SpectrumCount;
import zhanglab.inputconvertor.module.TargetDecoyAnalysis;
import zhanglab.inputconvertor.module.ToAutoRTInput;
import zhanglab.inputconvertor.module.ToFeatures;
import zhanglab.inputconvertor.module.ToMS2PIPInput;
import zhanglab.inputconvertor.module.ToTranslator;
import zhanglab.inputconvertor.module.NetMHCpan;
import zhanglab.inputconvertor.module.TopXgInput;

public class Run {

	public static void main(String[] args) throws IOException, ParseException {
		long startTime = System.currentTimeMillis();
		System.out.println("InputConvetor v0.0.0");
		
		Options options = new Options();
		
		// Options
		options.addOption("m", "module", true, "mode (ex> pxg_input, ms2pip_input, autort_input");
		options.addOption("t", "tool", true, "input tool (ex> casanovo, pnovo3, peaks, pxg)");
		
		options.addOption("i", "input", true, "A path of folder including .csnv.mztab/.res files");
		options.addOption("f", "mgf", true, "A path of folder including .mgf files");
		options.addOption("p", "pattern", true, "Batch pattern. It assumes that PSMs with a title containing the batch pattern will be recognized as the same batch");
		options.addOption("o", "output_prefix", true, "Output prefix");
		
		options.addOption("C", "charge", true, "ic_charge or tool-reported charge index. Use ic_charge for pNovo3, otherwise use tool-reported charge.");
		options.addOption("S", "score", true, "score index");
		
		options.addOption("M", "ms2pip", true, "MS2PIP output mgf.");
		options.addOption("A", "autort", true, "AutoRT output file.");
		options.addOption("N", "netmhcpan", true, "NetMHCpan4.1 output file.");
		options.addOption("T", "frag_tol", true, "Fragment tolerance (Da).");
		options.addOption("E", "ppm_err", true, "PPM error index. If not specified, then automatically calculates the PPM error.");
		
		options.addOption("d", "fdr", true, "FDR ratio. ex> 0.01, 0.05...");
		
		CommandLineParser parser = new DefaultParser();
        CommandLine cmd = parser.parse(options, args);
        
        // parse input
        String mode = cmd.getOptionValue("m");
        String inputTool = cmd.getOptionValue("t");

        
        if(mode.equalsIgnoreCase("pxg_input")) {
    		// args: -i -p -b -f
        	TopXgInput topXgInput = null;
        	if(inputTool.equalsIgnoreCase(InputConvertorConstants.CASANOVO)) {
        		topXgInput = new Casanovo();
        	} else if(inputTool.equalsIgnoreCase(InputConvertorConstants.PNOVO3)) {
        		topXgInput = new pNovo3();
        	} else if(inputTool.equalsIgnoreCase(InputConvertorConstants.PEAKS)) {
        		topXgInput = new PEAKS();
        	}
        	
        	topXgInput.topXgInputFormat(cmd);
        }
        
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
        
        else if(mode.equalsIgnoreCase("autort_input")) {
        	// args: -i -p -S
        	// -i ./
        	// -p C3N-01488.T.csnv
        	// --score 11
        	ToAutoRTInput toAutoRTInput = null;
        	if(inputTool.equalsIgnoreCase(InputConvertorConstants.PXG)) {
        		toAutoRTInput = new pXg();
        	} else if(inputTool.equalsIgnoreCase(InputConvertorConstants.COMET)) {
        		toAutoRTInput = new Comet();
        	}
        	
        	toAutoRTInput.toAutoRTInputFormat(cmd);
        	
        }
        else if(mode.equalsIgnoreCase("feature_gen")) {
        	ToFeatures toFeatures = null;
        	// -i -p -f -M -A -E -T -C -S
        	toFeatures = new pXg();
        	
        	toFeatures.toFeatures(cmd);
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
        	
        	// -m -t -i -p -o
        	if(inputTool.equalsIgnoreCase(InputConvertorConstants.ARRIBA)) {
        		toTranslator = new Arriba();
        	}
        	
        	toTranslator.doTranslate(cmd);
        }
        
        long endTime = System.currentTimeMillis();
        
        System.out.println((endTime-startTime)/1000+" sec");
		
	}
}
