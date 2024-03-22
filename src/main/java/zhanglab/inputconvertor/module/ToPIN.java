package zhanglab.inputconvertor.module;

import java.io.BufferedReader;
import java.io.BufferedWriter;
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

import zhanglab.inputconvertor.data.Peptide;
import zhanglab.inputconvertor.env.InputConvertorConstants;

public class ToPIN {
	
	public static String inputFilePath;
	public static String outputFilePath;
	public static String tool;
	
	private ToPIN() {};
	
	public static void convertToPIN (String[] args) throws IOException {
		parseOptions(args);
		
		if(tool.equalsIgnoreCase(InputConvertorConstants.COMET)) {
			parseComet();
		} else if(tool.equalsIgnoreCase(InputConvertorConstants.PXG)) {
			parsepXg();
		}
	}
	
	public static void parsepXg () throws IOException {
		BufferedReader BR = new BufferedReader(new FileReader(inputFilePath));
		BufferedWriter BW = new BufferedWriter(new FileWriter(outputFilePath));
		
		String[] pXgHeaders = BR.readLine().split("\t");
		
		int deltaScoreIdx = -1,
		log2ReadsIdx = -1, log2MQIdx = -1, genomicIDIdx = -1,
		ppmIdx = -1,
		saIdx = -1,
		deltaRTIdx = -1,
		specIdIdx = -1,
		labelIdx = -1, scoreIdx = -1, chargeIdx = -1;
		
		int icPeptideIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_PEPTIDE_FIELD_NAME);
		int infPeptideIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_INFERRED_PEPTIDE_FIELD_NAME);
		
		specIdIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_SPEC_ID_FEILD_NAME);
		genomicIDIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_GENOMIC_ID_FEILD_NAME);
		labelIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_LABEL_FEILD_NAME);
		scoreIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_SEARCH_SCORE_FIELD_NAME);
		log2ReadsIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_READ_FIELD_NAME);
		log2MQIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_MQ_FIELD_NAME);
		ppmIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_PPM_FIELD_NAME);
		saIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_SA_FIELD_NAME);
		deltaRTIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_DELTA_RT_FIELD_NAME);
		deltaScoreIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_DELTA_SCORE_FIELD_NAME);
		chargeIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_CHARGE_FIELD_NAME);
		
		int maxCharge = 0;
		int minCharge = 100;
		Hashtable<String, Integer> scanNums = new Hashtable<String, Integer>();
		ArrayList<String> records = new ArrayList<String>();
		String line = null;
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			int charge = (int) (Double.parseDouble(fields[chargeIdx]));
			
			maxCharge = Math.max(maxCharge, charge);
			minCharge = Math.min(minCharge, charge);
			
			records.add(line);
			
		}
		
		BR.close();
		
		// Build PIN File
		BW.append("SpecId").append("\t")
		.append("Label").append("\t")
		.append("ScanNr").append("\t")
		.append("MainScore").append("\t")
		.append("DeltaScore").append("\t")
		.append("Log2Reads").append("\t")
		.append("Log2MeanQScore").append("\t")
		.append("Length=7").append("\t")
		.append("Length=8").append("\t")
		.append("Length>8").append("\t")
		.append("abs(PPMError)").append("\t")
		.append("SA").append("\t")
		.append("abs(deltaRT)").append("\t");
		
		for(int i=minCharge; i<=maxCharge; i++) {
			BW.append("Charge"+i+"\t");
		}
		
		BW.append("Peptide").append("\t")
		.append("Proteins");
		
		BW.newLine();
		
		for(int i=0; i<records.size(); i++) {
			String[] fields = records.get(i).split("\t");
			String specId = fields[specIdIdx];
			
			Integer scanNum = scanNums.get(specId);
			if(scanNum == null) {
				scanNum = scanNums.size()+1;
				scanNums.put(specId, scanNum);
			}
			// skip if there is no SA/DeltaRT
			if(fields[saIdx].length() == 0 || fields[deltaRTIdx].length() ==0 ) {
				continue;
			}
			
			Peptide peptide = new Peptide(fields[icPeptideIdx], InputConvertorConstants.IC_CONVERTOR);
			peptide.icPeptideToILDetermination(fields[infPeptideIdx]);
			double log2Reads = Math.log(Double.parseDouble(fields[log2ReadsIdx]) +1) / Math.log(2);
			double log2MQ = Math.log(Double.parseDouble(fields[log2MQIdx])) / Math.log(2);
			int charge = (int) (Double.parseDouble(fields[chargeIdx]));
			String genomicID = fields[labelIdx].equalsIgnoreCase("-1") ? "XXX_"+fields[genomicIDIdx] : fields[genomicIDIdx];
			int length7 = peptide.stripPeptide.length() == 7 ? 1 : 0;
			int length8 = peptide.stripPeptide.length() == 8 ? 1 : 0;
			int lengthGreater8 = peptide.stripPeptide.length() > 8 ? 1 : 0;
			
			
			BW.append(fields[specIdIdx]).append("\t")
			.append(fields[labelIdx]).append("\t")
			.append(scanNum+"").append("\t")
			.append(fields[scoreIdx]).append("\t")
			.append(fields[deltaScoreIdx]).append("\t")
			.append(log2Reads+"").append("\t")
			.append(log2MQ+"").append("\t")
			.append(length7+"").append("\t")
			.append(length8+"").append("\t")
			.append(lengthGreater8+"").append("\t")
			.append(fields[ppmIdx]).append("\t")
			.append(fields[saIdx]).append("\t")
			.append(fields[deltaRTIdx]).append("\t");
			
			for(int c=minCharge; c<=maxCharge; c++) {
				if(c == charge) {
					BW.append("1").append("\t");
				} else {
					BW.append("0").append("\t");
				}
			}
			
			BW.append(peptide.modPeptide).append("\t")
			.append(genomicID);
			
			BW.newLine();
		}
		
		BW.close();
	}
	
	public static void parseComet () throws IOException {
		BufferedReader BR = new BufferedReader(new FileReader(inputFilePath));
		BufferedWriter BW = new BufferedWriter(new FileWriter(outputFilePath));
		
		String line = null;
		
		String pinHeader = BR.readLine().split("###")[1].trim();
		String[] pinHeaderSplit = pinHeader.split("\t");
		
		int peptideIdx = InputConvertorConstants.getFieldIndex(pinHeaderSplit, InputConvertorConstants.COMET_PIN_PEPTIDE_FIELD_NAME);
		int proteinIdx = InputConvertorConstants.getFieldIndex(pinHeaderSplit, InputConvertorConstants.COMET_PIN_PROTEIN_FIELD_NAME);
		int icPPMIdx = InputConvertorConstants.getFieldIndex(pinHeaderSplit, InputConvertorConstants.IC_PPM_FIELD_NAME);
		int icSAIdx = InputConvertorConstants.getFieldIndex(pinHeaderSplit, InputConvertorConstants.IC_SA_FIELD_NAME);
		int icDeltaRTIdx = InputConvertorConstants.getFieldIndex(pinHeaderSplit, InputConvertorConstants.IC_DELTA_RT_FIELD_NAME);
		
		// ban features
		int mass1Idx = InputConvertorConstants.getFieldIndex(pinHeaderSplit, InputConvertorConstants.COMET_PIN_MASS_FIELD_NAME);
		int mass2Idx = InputConvertorConstants.getFieldIndex(pinHeaderSplit, InputConvertorConstants.COMET_PIN_CALCMASS_FIELD_NAME);
		int mass3Idx = InputConvertorConstants.getFieldIndex(pinHeaderSplit, InputConvertorConstants.COMET_PIN_EXPMASS_FIELD_NAME);
		int mass4Idx = InputConvertorConstants.getFieldIndex(pinHeaderSplit, InputConvertorConstants.COMET_PIN_DM_FIELD_NAME);
		int mass5Idx = InputConvertorConstants.getFieldIndex(pinHeaderSplit, InputConvertorConstants.COMET_PIN_ABSDM_FIELD_NAME);
		int enzNIdx = InputConvertorConstants.getFieldIndex(pinHeaderSplit, InputConvertorConstants.COMET_PIN_ENZN_FIELD_NAME);
		int enzCIdx = InputConvertorConstants.getFieldIndex(pinHeaderSplit, InputConvertorConstants.COMET_PIN_ENZC_FIELD_NAME);
		int enzIntIdx = InputConvertorConstants.getFieldIndex(pinHeaderSplit, InputConvertorConstants.COMET_PIN_ENZ_INT_FIELD_NAME);
		
		int[] banFeatures = {mass1Idx, mass2Idx, mass3Idx, mass4Idx, mass5Idx,
							enzNIdx, enzCIdx, enzIntIdx};
		
		BW.append(pinHeaderSplit[0]);
		
		for(int i=1; i<peptideIdx; i++) {
			boolean isBan = false;
			for(int j=0; j<banFeatures.length; j++) {
				if(banFeatures[j] == i) {
					isBan = true;
				}
			}
			if(isBan) continue;
			BW.append("\t").append(pinHeaderSplit[i]);
		}
		BW.append("\t").append(pinHeaderSplit[icPPMIdx])
		.append("\t").append(pinHeaderSplit[icSAIdx])
		.append("\t").append(pinHeaderSplit[icDeltaRTIdx])
		.append("\t").append(pinHeaderSplit[peptideIdx])
		.append("\t").append(pinHeaderSplit[proteinIdx]);
		BW.newLine();
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("###")[1].trim().split("\t");
			
			BW.append(fields[0]);
			for(int i=1; i<peptideIdx; i++) {
				boolean isBan = false;
				for(int j=0; j<banFeatures.length; j++) {
					if(banFeatures[j] == i) {
						isBan = true;
					}
				}
				if(isBan) continue;
				
				BW.append("\t").append(fields[i]);
			}
			BW.append("\t").append(fields[icPPMIdx])
			.append("\t").append(fields[icSAIdx])
			.append("\t").append(fields[icDeltaRTIdx])
			.append("\t").append(fields[peptideIdx])
			.append("\t").append(fields[proteinIdx]);
			BW.newLine();
		}
		
		BR.close();
		BW.close();
	}
	
	public static void parseOptions (String[] args) {
		CommandLine cmd = null;
		Options options = new Options();
		
		// Mandatory
		Option optionInputTxt = Option.builder("i")
				.longOpt("input").argName("feat")
				.hasArg()
				.required(true)
				.desc("feat file")
				.build();
		
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("pin")
				.hasArg()
				.required(true)
				.desc("pin output file")
				.build();
		
		Option optionTool = Option.builder("t")
				.longOpt("tool").argName("comet/pxg")
				.hasArg()
				.required(false)
				.desc("specify a format of result to convert")
				.build();
		
		options.addOption(optionInputTxt)
		.addOption(optionTool)
		.addOption(optionOutput);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
	    ArrayList<String> tmpArgs = new ArrayList<String>();
	    for(int i=0; i<args.length; i++) {
	    	if(args[i].equalsIgnoreCase("-i") || args[i].equalsIgnoreCase("--input") ||
			args[i].equalsIgnoreCase("-o") || args[i].equalsIgnoreCase("--output") ||
			args[i].equalsIgnoreCase("-t") || args[i].equalsIgnoreCase("--tool")) {
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
