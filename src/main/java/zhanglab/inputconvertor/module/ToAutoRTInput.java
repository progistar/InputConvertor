package zhanglab.inputconvertor.module;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.data.AutoRTRecord;
import zhanglab.inputconvertor.data.Peptide;
import zhanglab.inputconvertor.env.InputConvertorConstants;

public class ToAutoRTInput {
	public static String inputFilePath = null;
	public static String outputFilePath = null;
	
	private ToAutoRTInput() {};
	
	public static void toAutoRTInputFormat (String[] args, String tool) throws IOException, ParseException {
		parseOptions(args);
		
		
        int scoreIdx		= -1;
        int icRTIdx			= -1;
        int icPeptideIdx 	= -1;
        int infPeptideIdx	= -1;
        
        File inputFile = new File(inputFilePath);
        BufferedReader BR = new BufferedReader(new FileReader(inputFile));
		System.out.println("Read "+inputFile.getName());
		String line = null;
		String pXgHeader = BR.readLine(); // skip header
		
		// get index
		String[] pXgHeaderSplit = pXgHeader.split("\t");
		scoreIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_SEARCH_SCORE_FIELD_NAME);
		icRTIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_RT_FIELD_NAME);
		icPeptideIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_PEPTIDE_FIELD_NAME);
		
		if(tool.equalsIgnoreCase(InputConvertorConstants.PXG)) {
			infPeptideIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_INFERRED_PEPTIDE_FIELD_NAME);
    	} else if(tool.equalsIgnoreCase(InputConvertorConstants.COMET)) {
    		// TODO
    		// parse for COMET!
    		System.exit(1);
    		infPeptideIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.COMET_PEPTIDE_FIELD_NAME);
    	}
		
		BufferedWriter BW = new BufferedWriter(new FileWriter(outputFilePath));
		// building AutoRT header ///////////////////////////////////////////////
		
		String header = InputConvertorConstants.AUTORT_HEADER_X+"\t"+InputConvertorConstants.AUTORT_HEADER_Y;
		BW.append(header);
		BW.newLine();
		
		/////////////////////////////////////////////////////////////////////////
		
		
		ArrayList<AutoRTRecord> records = new ArrayList<AutoRTRecord>();
		int discardedPRSMs = 0;
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			
			double score = Double.parseDouble(fields[scoreIdx]);
			double rt	 = Double.parseDouble(fields[icRTIdx]);
			Peptide icPeptide	=	new Peptide(fields[icPeptideIdx], InputConvertorConstants.IC_CONVERTOR);
			String infPeptide	=	fields[infPeptideIdx];
			
			AutoRTRecord record = new AutoRTRecord();
			
			record.score = score;
			// sec to min
			record.rt = (rt/60)+"";
			// IL determination from infPeptide
			icPeptide.icPeptideToILDetermination(infPeptide);
			
			/****************************************************************************
			 * TODO: variable modifications
			 * Current version ignores other PTMs except for M_OXI.
			 * Just drop the records!
			 ****************************************************************************/
			//check if there are other PTMs...
			String checkPeptide = icPeptide.stripPeptide.replaceAll("[+-0123456789.*]", "");
			if(!checkPeptide.equalsIgnoreCase(icPeptide.stripPeptide)) {
				// drop the peptide
				discardedPRSMs++;
				continue;
			}
			
			// convert IC_M_OXI to AUTORT_M_OXI
			// e.g. ACM+15.995GG => AC1GG
			icPeptide.toAutoRTModPeptide();
			
			record.modifiedPeptide = icPeptide.modPeptide;
			record.fullRecord = line;
			
			// add record
			records.add(record);
			record.idx = records.size();
		}
		
		System.out.println("A total of PRSMs: "+records.size());
		if(discardedPRSMs != 0) {
			System.out.println("Note! "+discardedPRSMs+" PRSMs were discarded due to having unsupported modifications");
		}
		
		// remove duplicated peptides and remain top-scored RT per peptide.
		Hashtable<String, AutoRTRecord> hasPeptide = new Hashtable<String, AutoRTRecord>();
		for(int i=0; i<records.size(); i++) {
			AutoRTRecord record = records.get(i);
			AutoRTRecord hasRecord = hasPeptide.get(record.modifiedPeptide);
			
			if(hasRecord == null) {
				hasPeptide.put(record.modifiedPeptide, record);
			} else {
				if(record.score > hasRecord.score) {
					hasPeptide.put(record.modifiedPeptide, record);
				}
			}
		}
		System.out.println("Remain top-scored PRSM per peptide ... ");
		//// clear the records and insult again from hasPeptide
		records.clear();
		hasPeptide.forEach((p, record)->{
			records.add(record);
		});
		// sort by index to make sure the order.
		Collections.sort(records);
		System.out.println("A total of peptides to be considered in AutoRT: "+records.size());
		
		for(int i=0; i<records.size(); i++) {
			AutoRTRecord record = records.get(i);
			BW.append(record.modifiedPeptide+"\t"+record.rt);
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
				.longOpt("input").argName("res")
				.hasArg()
				.required(true)
				.desc("casanovo output file")
				.build();
		
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("file path")
				.hasArg()
				.required(true)
				.desc("output file path")
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
		    cmd = parser.parse(options, nArgs, true);
		    
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
