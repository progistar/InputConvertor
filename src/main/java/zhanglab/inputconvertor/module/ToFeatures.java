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

import zhanglab.inputconvertor.data.DeepLCRecord;
import zhanglab.inputconvertor.data.MS2PIPRecord;
import zhanglab.inputconvertor.data.Peptide;
import zhanglab.inputconvertor.data.Spectra;
import zhanglab.inputconvertor.data.Spectrum;
import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.env.ProteomeConstants;
import zhanglab.inputconvertor.function.CalculateSA;

class FeatureRecord implements Comparable <FeatureRecord> {
	public String[] record;
	public double searchScore;
	public String specId;
	public String label;
	public double deltaRT = Double.MAX_VALUE;
	public double ppmError = Double.MAX_VALUE;
	
	@Override
	public int compareTo(FeatureRecord o) {
		if(this.searchScore < o.searchScore) {
			return 1;
		} else if(this.searchScore > o.searchScore) {
			return -1;
		}
		return 0;
	}
}

public class ToFeatures {
	public static String inputFilePath = null;
	public static String outputFilePath = null;
	public static String spectrumFilePath = null;
	public static double tol = .0;
	public static String deepLCFilePath = null;
	public static String predictedFilePath = null;
	
	
	private ToFeatures() {};
	
	public static void toFeatureFormat (String[] args) throws IOException, ParseException {
		parseOptions(args);
		
		
        int scoreIdx		= -1;
        int icPeptideIdx 	= -1;
        int infPeptideIdx	= -1;
        int chargeIdx	 		= -1;
        int specIdIdx		= -1;
        int labelIdx		= -1;
		
		// read deepLC output file  //////////////////////////////////////////////////
        BufferedReader BRdeepLC = new BufferedReader(new FileReader(deepLCFilePath));
        System.out.println("Read deepLC...: "+deepLCFilePath);
        ArrayList<DeepLCRecord> deepLCRecords = new ArrayList<DeepLCRecord>();
        Hashtable<String, Double> icPeptideToDeltaRT = new Hashtable<String, Double>();
        BRdeepLC.readLine(); // skip header
        String line = null;
        
        while((line = BRdeepLC.readLine()) != null) {
        	DeepLCRecord deepLCRecord = new DeepLCRecord();
        	deepLCRecords.add(deepLCRecord);
        	
        	String[] fields = line.split("\t");
        	
        	deepLCRecord.idx = deepLCRecords.size();
        	deepLCRecord.fullRecord = line;
        	deepLCRecord.modifiedPeptide = fields[InputConvertorConstants.DEEPLC_PEPTIDE_IDX];
        	deepLCRecord.rt = fields[InputConvertorConstants.DEEPLC_RT_IDX];
        	deepLCRecord.predRT = fields[InputConvertorConstants.DEEPLC_PRED_IDX];
        	
        	Peptide peptide = new Peptide(deepLCRecord.modifiedPeptide, InputConvertorConstants.PROSIT);
        	
        	// calculate |delta RT|
        	double deltaRT = Double.parseDouble(deepLCRecord.predRT) - Double.parseDouble(deepLCRecord.rt);
        	// deltaRT = Math.abs(deltaRT);
        	icPeptideToDeltaRT.put(peptide.modPeptide, deltaRT);
        }
        System.out.println("The number of deepLC peptides: "+deepLCRecords.size());
        System.out.println("The number of hashed peptides: "+icPeptideToDeltaRT.size());
        BRdeepLC.close();
        ///////////////////// ///////////////// End of reading deepLC output file ////

		File inputFile = new File(inputFilePath);
        
		BufferedReader BR = new BufferedReader(new FileReader(inputFile));
		System.out.println("read: "+inputFile.getName());
		String header = BR.readLine(); // read header
		header += "\t" + InputConvertorConstants.IC_PPM_FIELD_NAME + "\t" + InputConvertorConstants.IC_SA_FIELD_NAME + "\t" + InputConvertorConstants.IC_DELTA_RT_FIELD_NAME;
		header += "\tDeltaRT(original)\tDeltaRT(shifted)\tppmError(original)\tppmError(shifted)";
		// building header ///////////////////////////////////////////
    	// get index
		String[] pXgHeaderSplit = header.split("\t");
		pXgHeaderSplit = header.split("\t");
		
		scoreIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_SEARCH_SCORE_FIELD_NAME);
		icPeptideIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_PEPTIDE_FIELD_NAME);
		infPeptideIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_INFERRED_PEPTIDE_FIELD_NAME);
		chargeIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_CHARGE_FIELD_NAME);
		specIdIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_SPEC_ID_FEILD_NAME);
		labelIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_LABEL_FEILD_NAME);
		
		int ppmIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_PPM_FIELD_NAME);
		int saIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_SA_FIELD_NAME);
		int deltaRTIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_DELTA_RT_FIELD_NAME);
        ////////////////////////////////// End of building header ////////////////
		
        String predictType = InputConvertorConstants.PROSIT;
        Spectra predictedSpectra = new Spectra(predictedFilePath, predictType);
        Spectra spectra = new Spectra(spectrumFilePath, InputConvertorConstants.EXPERIMENTAL_SPECTRUM);
        
        ArrayList<FeatureRecord> newRecords = new ArrayList<FeatureRecord>();
        
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String[] newFields = new String[pXgHeaderSplit.length];
			
			for(int i=0; i<newFields.length; i++) {
				if(fields.length <= i) {
					newFields[i] = "";
				} else {
					newFields[i] = fields[i];
				}
			}
			
			// Building record ////////////////////////////////////////
			String specID = fields[specIdIdx];
			String title = specID.split("\\|")[0];
			
			
			FeatureRecord newRecord = new FeatureRecord();
			newRecord.record = newFields;
			newRecord.searchScore = Double.parseDouble(newFields[scoreIdx]);
			newRecord.specId = specID;
			newRecord.label = newFields[labelIdx];
			
			if(!spectra.isComplete) {
				System.out.println(spectrumFilePath);
				System.out.println("Sever error was occurred when reading spectrum!");
				System.exit(1);
			}
			
			String infPeptide = fields[infPeptideIdx];
			Peptide peptide = new Peptide(fields[icPeptideIdx], InputConvertorConstants.IC_CONVERTOR);
			peptide.icPeptideToILDetermination(infPeptide);
			String charge	= fields[chargeIdx];
			double chargeNum = Double.parseDouble(charge);
			
			String modifications = MS2PIPRecord.getModifications(peptide);
			String predictKey = MS2PIPRecord.getMS2PIPKey(peptide.stripPeptide, modifications, charge);
			
			
			Spectrum expSpectrum = spectra.scanIndexer.get(title);
			Spectrum predictedSpectrum = predictedSpectra.scanIndexer.get(predictKey);
			
			if(expSpectrum != null && predictedSpectrum != null) {
				// set peptide
				expSpectrum.peptide = peptide;
				predictedSpectrum.peptide = peptide;
				
				// MS2PIP dependent parameter here!! //////////////////////////////////
				// #TODO #WARNING /////////////////////////////////////////////////////
				// The current non-tryptic model in MS2PIP predicts only singly-charged
				// fragment ions.
				int chargeLimit = 0;
				if(predictType.equalsIgnoreCase(InputConvertorConstants.MS2PIP)) {
					chargeLimit = 1;
				} else if(predictType.equalsIgnoreCase(InputConvertorConstants.PROSIT)) {
					chargeLimit = 2;
				}
				
				double sa = CalculateSA.calculate(expSpectrum, predictedSpectrum, tol, (int) chargeNum,  chargeLimit);
				///////////////////////////////////////////////////////////////////////
				
				// Calculate delta Mz
				// isotope check
				double deltaMz = 100000;
				for(int iso=-1; iso<3; iso++) {
					double expMz = expSpectrum.precursorMz;
					double thrMz = (peptide.getMass(true)+chargeNum*ProteomeConstants.Proton + iso * ProteomeConstants.IsotopeSpace) / chargeNum;
					double thisDelta = ((expMz - thrMz) / thrMz) * Math.pow(10, 6);
					
					if(Math.abs(thisDelta) < Math.abs(deltaMz)) {
						deltaMz = thisDelta;
					}
				}
				
				
				newFields[ppmIdx] = "" + deltaMz;
				newFields[saIdx] = "" + sa;
				
				newRecord.ppmError = deltaMz;
				

				// Calculate delta RT
				if(icPeptideToDeltaRT.get(peptide.modPeptide) != null) {
					double deltaRT = icPeptideToDeltaRT.get(peptide.modPeptide);
					newFields[deltaRTIdx] = "" + deltaRT;
					newRecord.deltaRT = deltaRT;
				}
				
				newRecords.add(newRecord);
			}
			
			
			/////////////////////////////////// End of building record ///////////////////
		}
		BR.close();
		
		///////////// Top 100 median shift per spectrum file /////////
		topMedianShift(newRecords, deltaRTIdx, ppmIdx);
		
		File outputFile = new File(outputFilePath);
		boolean isAppended = outputFile.exists();
		BufferedWriter BW = new BufferedWriter(new FileWriter(outputFile, isAppended));
		
		if(!isAppended) {
			BW.append(header);
			BW.newLine();
		}
		for(FeatureRecord record :newRecords) {
			String[] fields = record.record;
			for(int i=0; i<fields.length; i++) {
				if(i != 0) {
					BW.append("\t");
				}
				BW.append(fields[i]);
			}
			BW.newLine();
		}
		BW.close();
	}
	
	/**
	 * Select target PSMs up to top 100.
	 * Calculate their median, and apply all shift by the median value.
	 * 
	 * @param records
	 * @param deltaRTIdx
	 */
	public static void topMedianShift (ArrayList<FeatureRecord> records, int deltaRTIdx, int ppmErrorIdx) {
		Collections.sort(records);
		Hashtable<String, String> isOverlappedSpectrum = new Hashtable<String, String>();
		ArrayList<Double> deltaRTs = new ArrayList<Double>();
		ArrayList<Double> ppmErrors = new ArrayList<Double>();
		
		for(int i=0; i<records.size(); i++) {
			FeatureRecord record = records.get(i);
			if(record.label.equalsIgnoreCase("-1")) {
				continue;
			}
			
			if(isOverlappedSpectrum.get(record.specId) == null) {
				isOverlappedSpectrum.put(record.specId, "");
				deltaRTs.add(record.deltaRT);
				ppmErrors.add(record.ppmError);
			}
			
			if(deltaRTs.size() == 1000) {
				break;
			}
		}
		int size = deltaRTs.size();
		double medianRT = deltaRTs.get(size/2) + deltaRTs.get((size-1)/2);
		double medianPPM = ppmErrors.get(size/2) + ppmErrors.get((size-1)/2);
		medianRT /= 2;
		medianPPM /= 2;
		
		double meanRT = 0;
		double meanPPM = 0;
		for(int i=0; i<size; i++) {
			meanRT += deltaRTs.get(i);
			meanPPM += ppmErrors.get(i);
		}
		meanRT /= size;
		meanPPM /= size;
		
		double stdRT = 0;
		double stdPPM = 0;
		for(int i=0; i<size; i++) {
			stdRT += Math.pow(deltaRTs.get(i) - meanRT,2);
			stdPPM += Math.pow(ppmErrors.get(i) - meanPPM,2);
		}
		stdRT /= size;
		stdPPM /= size;
		stdRT = Math.sqrt(stdRT);
		stdPPM = Math.sqrt(stdPPM);
		
		for(int i=0; i<records.size(); i++) {
			FeatureRecord record = records.get(i);
			
			record.record[deltaRTIdx+1] = record.deltaRT+"";
			record.record[deltaRTIdx+2] = (record.deltaRT - meanRT)/stdRT+"";
			
			record.record[deltaRTIdx+3] = record.ppmError+"";
			record.record[deltaRTIdx+4] = (record.ppmError - meanPPM)/stdPPM+"";
			
			record.deltaRT = Math.abs((record.deltaRT - meanRT)/stdRT);
			record.ppmError = Math.abs((record.ppmError - meanPPM)/stdPPM);
			record.record[deltaRTIdx] = record.deltaRT+"";
			record.record[ppmErrorIdx] = record.ppmError+"";
		}
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
		
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("pin")
				.hasArg()
				.required(true)
				.desc("pin file")
				.build();
		
		Option optionSpectrumFile = Option.builder("s")
				.longOpt("spectrum").argName("mgf/mzml")
				.hasArg()
				.required(true)
				.desc("spectrum file")
				.build();
		
		Option optionDeepLCFile = Option.builder("d")
				.longOpt("deeplc").argName("tsv")
				.hasArg()
				.required(true)
				.desc("deepLC output file")
				.build();
		
		Option optionPredictedFileP = Option.builder("p")
				.longOpt("pspectrum").argName("mgf")
				.hasArg()
				.required(true)
				.desc("predicted spectrum file")
				.build();
		
		Option optionTolerance = Option.builder("t")
				.longOpt("tolerance").argName("float")
				.hasArg()
				.required(true)
				.desc("fragment tolerance")
				.build();
		
		options.addOption(optionInput)
		.addOption(optionSpectrumFile)
		.addOption(optionDeepLCFile)
		.addOption(optionPredictedFileP)
		.addOption(optionTolerance)
		.addOption(optionOutput);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
	    ArrayList<String> tmpArgs = new ArrayList<String>();
	    for(int i=0; i<args.length; i++) {
	    	if(args[i].equalsIgnoreCase("-i") || args[i].equalsIgnoreCase("--input") ||
			args[i].equalsIgnoreCase("-s") || args[i].equalsIgnoreCase("--spectrum") ||
			args[i].equalsIgnoreCase("-d") || args[i].equalsIgnoreCase("--deeplc") ||
			args[i].equalsIgnoreCase("-p") || args[i].equalsIgnoreCase("--pspectrum") ||
			args[i].equalsIgnoreCase("-t") || args[i].equalsIgnoreCase("--tolerance") ||
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
		    
		    if(cmd.hasOption("d")) {
		    	deepLCFilePath = cmd.getOptionValue("d");
		    }
		    
		    if(cmd.hasOption("p")) {
		    	predictedFilePath = cmd.getOptionValue("p");
		    }
		    
		    if(cmd.hasOption("t")) {
		    	tol = Double.parseDouble(cmd.getOptionValue("t"));
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
