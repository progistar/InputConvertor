package zhanglab.inputconvertor.input;

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
import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.data.AutoRTRecord;
import zhanglab.inputconvertor.data.MS2PIPRecord;
import zhanglab.inputconvertor.data.NetMHCpanResult;
import zhanglab.inputconvertor.data.Peptide;
import zhanglab.inputconvertor.data.Spectra;
import zhanglab.inputconvertor.data.Spectrum;
import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.function.CalculateSA;
import zhanglab.inputconvertor.function.NetMHCpanParser;
import zhanglab.inputconvertor.module.ToAutoRTInput;
import zhanglab.inputconvertor.module.ToFeatures;
import zhanglab.inputconvertor.module.ToMS2PIPInput;
import zhanglab.inputconvertor.module.NetMHCpan;


public class pXg implements ToAutoRTInput, ToMS2PIPInput, ToFeatures, NetMHCpan {
	public pXg () {}
	
	@Override
	public void addNetMHCpanOutput(CommandLine cmd) throws IOException, ParseException {
		String inputFile = cmd.getOptionValue("i");
        String inputPattern = cmd.getOptionValue("p");

        int infPeptideIdx		= -1;
        
        File[] files = new File(inputFile).listFiles();
        if(files == null) {
        	System.out.println("-i option should be a path of folder including .fdr files.");
        	System.exit(1);
        }
        
        for(File file : files) {
        	if(file.getName().startsWith(".")) continue;
        	if(file.getName().endsWith(".fdr") && file.getName().contains(inputPattern)) {
        		
        		

                // Read NetMHCpan Output ////////////////////////////////////////////////////////
                NetMHCpanResult result = NetMHCpanParser.parseNetMHCpan(file.getAbsolutePath().replace(".fdr", ".netMHCpan.xls"));
                /////////////////////////////////// End of NetMHCpan Output reader //////////////
                
        		
        		String outputName = inputFile+"/"+file.getName().replace(".fdr", ".fdr.ba");
        		BufferedWriter BW = new BufferedWriter(new FileWriter(outputName));
        		BufferedReader BR = new BufferedReader(new FileReader(file));
        		String line = null;
        		String header = BR.readLine();
        		String[] headerSplit = header.split("\t");
        		
        		infPeptideIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.IC_INFERREND_PEPTIDE_FIELD_NAME);
        		
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
        }
	}

	@Override
	public void toNetMHCpanInputFormat(CommandLine cmd) throws IOException, ParseException {
		String inputFile = cmd.getOptionValue("i");
        String inputPattern = cmd.getOptionValue("p");

        int infPeptideIdx		= -1;
        
        File[] files = new File(inputFile).listFiles();
        if(files == null) {
        	System.out.println("-i option should be a path of folder including .fdr files.");
        	System.exit(1);
        }
        
        for(File file : files) {
        	if(file.getName().startsWith(".")) continue;
        	if(file.getName().endsWith(".fdr") && file.getName().contains(inputPattern)) {
        		
        		String outputName = inputFile+"/"+file.getName().replace(".fdr", ".netMHCpan.input");
        		BufferedWriter BW = new BufferedWriter(new FileWriter(outputName));
        		BufferedReader BR = new BufferedReader(new FileReader(file));
        		
        		String line = null;
        		String[] header = BR.readLine().split("\t");
        		
        		infPeptideIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.IC_INFERREND_PEPTIDE_FIELD_NAME);
        		
        		Hashtable<String, String> isDuplicated = new Hashtable<String, String>();
        		while((line = BR.readLine()) != null) {
        			String[] fields = line.split("\t");
        			String infPeptide = fields[infPeptideIdx];
        			
        			if(infPeptide.length() >= 8 && infPeptide.length() <= 15) {

            			if(isDuplicated.get(infPeptide) == null) {
            				BW.append(infPeptide);
            				BW.newLine();
            				isDuplicated.put(infPeptide, "");
            			}
        			}
        			
        		}
        		
        		BR.close();
        		BW.close();
        		
        	}
        }
        
	}
	
	public void toFeatures (CommandLine cmd) throws IOException, ParseException {
		
		/////////////////////
        String inputFile = cmd.getOptionValue("i");
        String inputPattern = cmd.getOptionValue("p");
        String mgfFileBase = cmd.getOptionValue("f");
        String ms2pipOutputFile = cmd.getOptionValue("M");
        String autortOutputFile = cmd.getOptionValue("A");
        double tol = Double.parseDouble(cmd.getOptionValue("T"));
        int scoreIdx			= Integer.parseInt(cmd.getOptionValue("S"));
        int chargeIdx	 		= Integer.parseInt(cmd.getOptionValue("C"));
        
        int icPeptideIdx 		= -1;
        int infPeptideIdx		= -1;
        int ppmErrorIdx			= -1;
        
        if( cmd.getOptionValue("E") != null) {
        	ppmErrorIdx			= Integer.parseInt(cmd.getOptionValue("E"));
        }
        
        File[] files = new File(inputFile).listFiles();
        if(files == null) {
        	System.out.println("-i option should be a path of folder including .pXg files");
        	System.exit(1);
        }
        
        Spectra ms2pipSpectra = new Spectra(ms2pipOutputFile);
        
        // read AutoRT output file  //////////////////////////////////////////////////
        BufferedReader BRautoRT = new BufferedReader(new FileReader(autortOutputFile));
        ArrayList<AutoRTRecord> autoRTRecords = new ArrayList<AutoRTRecord>();
        Hashtable<String, Double> icPeptideToDeltaRT = new Hashtable<String, Double>();
        String line = BRautoRT.readLine(); // skip header
        
        while((line = BRautoRT.readLine()) != null) {
        	AutoRTRecord autoRTRecord = new AutoRTRecord();
        	autoRTRecords.add(autoRTRecord);
        	
        	String[] fields = line.split("\t");
        	
        	autoRTRecord.idx = autoRTRecords.size();
        	autoRTRecord.fullRecord = line;
        	autoRTRecord.modifiedPeptide = fields[InputConvertorConstants.AUTORT_X_IDX];
        	autoRTRecord.rt = fields[InputConvertorConstants.AUTORT_Y_IDX];
        	autoRTRecord.predRT = fields[InputConvertorConstants.AUTORT_PRED_Y_IDX];
        	
        	Peptide peptide = new Peptide(autoRTRecord.modifiedPeptide, InputConvertorConstants.AUTORT);
        	
        	// calculate |delta RT|
        	double deltaRT = Double.parseDouble(autoRTRecord.predRT) - Double.parseDouble(autoRTRecord.rt);
        	deltaRT = Math.abs(deltaRT);
        	
        	icPeptideToDeltaRT.put(peptide.modPeptide, deltaRT);
        }
        System.out.println("The number of AutoRT peptides: "+autoRTRecords.size());
        System.out.println("The number of hashed peptides: "+icPeptideToDeltaRT.size());
        BRautoRT.close();
        ///////////////////// ///////////////// End of reading AutoRT output file ////
        
        String header = null;
        for(File file : files) {
        	if(file.getName().startsWith(".")) continue;
        	if(file.getName().endsWith(".pXg") && file.getName().contains(inputPattern)) {
        		int skipRecords = 0;
        		
        		String outputName = inputFile+"/"+file.getName().replace(".pXg", ".pXg.feat");
        		BufferedWriter BW = new BufferedWriter(new FileWriter(outputName));
        		
        		System.out.println("read: "+file.getName());
        		BufferedReader BR = new BufferedReader(new FileReader(file));
        		line = BR.readLine(); // read header
        		
        		// building header ///////////////////////////////////////////
                if(header == null) {
                	header = line;
                	
                	// get index
            		String[] pXgHeaderSplit = header.split("\t");
            		icPeptideIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_PEPTIDE_FIELD_NAME);
            		infPeptideIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_INFERREND_PEPTIDE_FIELD_NAME);
                	
                	
                	BW.append(header).append("\t")
                	.append(InputConvertorConstants.IC_PPM_FIELD_NAME).append("\t")
    				.append(InputConvertorConstants.IC_SA_FIELD_NAME).append("\t")
    				.append(InputConvertorConstants.IC_DELTA_RT_FIELD_NAME);
    				BW.newLine();
                }
                ////////////////////////////////// End of building header ////////////////
        		
                Hashtable<String, Spectra> spectraHash = new Hashtable<String, Spectra>();
                int maxCharge = 0;
                int minCharge = 100;
                
        		while((line = BR.readLine()) != null) {
        			String[] fields = line.split("\t");
        			
        			
        			// Building record ////////////////////////////////////////
        			String specID = fields[0];
        			String title = specID.split("\\|")[0];
        			String fileName = title.split("\\.")[0];
        			
    				String mgfFile = mgfFileBase+"/"+fileName+".mgf";
    				
    				Spectra spectra = spectraHash.get(mgfFile);
    				if(spectra == null) {
    					spectra = new Spectra(mgfFile);
    					spectraHash.put(mgfFile, spectra);
    				}
    				
    				if(!spectra.isComplete) {
    					System.out.println(mgfFile);
    					System.out.println("Sever error was occurred when reading spectrum!");
    					System.exit(1);
    				}
    				
    				String infPeptide = fields[infPeptideIdx];
    				Peptide peptide = new Peptide(fields[icPeptideIdx], InputConvertorConstants.IC_CONVERTOR);
    				peptide.icPeptideToILDetermination(infPeptide);
    				String charge	= fields[chargeIdx];
    				double chargeNum = Double.parseDouble(charge);
    				
    				maxCharge = Math.max(maxCharge, (int)chargeNum);
    				minCharge = Math.min(maxCharge, (int)chargeNum);
    				
    				String modifications = MS2PIPRecord.getModifications(peptide);
    				String ms2pipKey = MS2PIPRecord.getMS2PIPKey(peptide.stripPeptide, modifications, charge);
    				
    				Spectrum expSpectrum = spectra.scanIndexer.get(title);
    				Spectrum ms2pipSpectrum = ms2pipSpectra.scanIndexer.get(ms2pipKey);
    				
    				if(ms2pipSpectrum == null) {
    					skipRecords++;
    					continue;
    				}
    				
    				// set peptide
    				expSpectrum.peptide = peptide;
    				ms2pipSpectrum.peptide = peptide;
    				
    				// MS2PIP dependent parameter here!! //////////////////////////////////
    				// #TODO #WARNING /////////////////////////////////////////////////////
    				// The current non-tryptic model in MS2PIP predicts only singly-charged
    				// fragment ions.
    				double sa = CalculateSA.calculate(expSpectrum, ms2pipSpectrum, tol, (int) chargeNum,  2);
    				///////////////////////////////////////////////////////////////////////
    				
    				// Calculate delta Mz
    				double absDeltaMz = -1;
    				if(ppmErrorIdx == -1) {
    					double expMz = expSpectrum.precursorMz;
        				double thrMz = ms2pipSpectrum.precursorMz;
        				absDeltaMz = ((expMz - thrMz) / thrMz) * Math.pow(10, 6);
    				} else {
    					absDeltaMz = Double.parseDouble(fields[ppmErrorIdx]);
    				}
    				
    				absDeltaMz = Math.abs(absDeltaMz);
    				
    				// Calculate delta RT
    				double absDeltaRT = icPeptideToDeltaRT.get(peptide.modPeptide);
    				
    				BW.append(line).append("\t")
    				.append(absDeltaMz+"").append("\t")
    				.append(sa+"").append("\t")
    				.append(absDeltaRT+"");
    				
    				BW.newLine();
    				/////////////////////////////////// End of building record ///////////////////
        		}
        		BR.close();
        		BW.close();
        		
        		System.out.println("Unspported spectra: "+skipRecords);
        		
        		// Build PIN File
        		BW = new BufferedWriter(new FileWriter(inputFile+"/"+file.getName().replace(".pXg", ".pXg.feat.pin")));
        		BR = new BufferedReader(new FileReader(outputName));
        		
        		BW.append("SpecId").append("\t")
        		.append("Label").append("\t")
        		.append("ScanNr").append("\t")
        		.append("MainScore").append("\t")
        		.append("DeltaScore").append("\t")
        		.append("Log2Reads").append("\t")
        		.append("Log2MeanQScore").append("\t")
        		.append("abs(PPMError)").append("\t")
        		.append("SA").append("\t")
        		.append("abs(deltaRT)").append("\t");
        		
        		for(int i=minCharge; i<=maxCharge; i++) {
        			BW.append("Charge"+i+"\t");
        		}
        		
        		BW.append("Peptide").append("\t")
        		.append("Proteins");
        		
        		BW.newLine();
        		
        		String[] pXgHeaders = BR.readLine().split("\t");
        		
        		int specIdIdx, labelIdx, mainScoreIdx, deltaScoreIdx = -1,
        		log2ReadsIdx = -1, log2MQIdx = -1, ppmIdx = -1, saIdx = -1, deltaRTIdx = -1, genomicIDIdx = -1;
        		
        		
        		icPeptideIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_PEPTIDE_FIELD_NAME);
        		infPeptideIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_INFERREND_PEPTIDE_FIELD_NAME);
        		
        		specIdIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_SPEC_ID_FEILD_NAME);
        		genomicIDIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_GENOMIC_ID_FEILD_NAME);
        		labelIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_LABEL_FEILD_NAME);
        		mainScoreIdx = scoreIdx;
        		log2ReadsIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_READ_FIELD_NAME);
        		log2MQIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_MQ_FIELD_NAME);
        		ppmIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_PPM_FIELD_NAME);
        		saIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_SA_FIELD_NAME);
        		deltaRTIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_DELTA_RT_FIELD_NAME);
        		deltaScoreIdx = InputConvertorConstants.getFieldIndex(pXgHeaders, InputConvertorConstants.IC_DELTA_SCORE_FIELD_NAME);
        		//chargeIdx ##
        		
        		Hashtable<String, Integer> scanNums = new Hashtable<String, Integer>();
        		while((line = BR.readLine()) != null) {
        			String[] fields = line.split("\t");
        			String specId = fields[specIdIdx];
        			
        			Integer scanNum = scanNums.get(specId);
        			if(scanNum == null) {
        				scanNum = scanNums.size()+1;
        				scanNums.put(specId, scanNum);
        			}
        			
        			
        			
        			Peptide peptide = new Peptide(fields[icPeptideIdx], InputConvertorConstants.IC_CONVERTOR);
        			peptide.icPeptideToILDetermination(fields[infPeptideIdx]);
        			double log2Reads = Math.log(Double.parseDouble(fields[log2ReadsIdx]) +1) / Math.log(2);
        			double log2MQ = Math.log(Double.parseDouble(fields[log2MQIdx])) / Math.log(2);
        			int charge = (int) (Double.parseDouble(fields[chargeIdx]));
        			String genomicID = fields[labelIdx].equalsIgnoreCase("-1") ? "XXX_"+fields[genomicIDIdx] : fields[genomicIDIdx];
        			
        			BW.append(fields[specIdIdx]).append("\t")
        			.append(fields[labelIdx]).append("\t")
        			.append(scanNum+"").append("\t")
        			
        			.append(fields[mainScoreIdx]).append("\t")
        			.append(fields[deltaScoreIdx]).append("\t")
        			.append(log2Reads+"").append("\t")
        			.append(log2MQ+"").append("\t")
        			.append(fields[ppmIdx]).append("\t")
        			.append(fields[saIdx]).append("\t")
        			.append(fields[deltaRTIdx]).append("\t");
        			
        			for(int i=minCharge; i<=maxCharge; i++) {
        				if(i == charge) {
        					BW.append("1").append("\t");
        				} else {
        					BW.append("0").append("\t");
        				}
        			}
        			
        			BW.append(peptide.modPeptide).append("\t")
        			.append(genomicID);
        			
        			BW.newLine();
        		}
        		
        		BR.close();
        		BW.close();
        	}
        }
        
        
	}
	
	
	public void toMS2PIPInputFormat (CommandLine cmd) throws IOException, ParseException {
		boolean isDebugMode = true;
		
		// Options //////////
		
		/**
		 * -i -p -P -I -C
		 * 
		 * -i pXg/
		 * -p Set01
		 * -P 12
		 * -I 20
		 * -C 24
		 * 
		 * For example>
		 * 
		 * target files:
		 * pXg/abc_Set01.pXg
		 * pXg/dfg_Set01.pXg
		 * pXg/ggg_set01.pXg // case sensitive
		 * ...
		 * 
		 * result files
		 * pXg/abc_Set01.ms2pip.input
		 * pXg/dfg_Set01.ms2pip.input
		 * 
		 */
		/////////////////////
		
        String inputFile	= cmd.getOptionValue("i");
        String filePattern	= cmd.getOptionValue("p");
        int chargeIdx		= Integer.parseInt(cmd.getOptionValue("C"));
        int icPeptideIdx 	= -1;
        int infPeptideIdx	= -1;
        
        
        File[] files = new File(inputFile).listFiles();
        if(files == null) {
        	System.out.println("-i option should be a path of folder including .pXg files");
        	System.exit(1);
        }
        
        for(File file : files) {
        	if(file.getName().startsWith(".")) continue;
        	if(file.getName().endsWith(".pXg") && file.getName().contains(filePattern)) {
        		
        		BufferedReader BR = new BufferedReader(new FileReader(file));
        		System.out.println("Read "+file.getName());
        		String line = null;
        		String pXgHeader = BR.readLine(); // skip header
        		
        		// get index
        		String[] pXgHeaderSplit = pXgHeader.split("\t");
        		icPeptideIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_PEPTIDE_FIELD_NAME);
        		infPeptideIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_INFERREND_PEPTIDE_FIELD_NAME);
        		
        		
        		String outputPath = file.getAbsolutePath().replace(".pXg", ".ms2pip.input");
        		BufferedWriter BW = new BufferedWriter(new FileWriter(outputPath));
        		// building MS2PIP header ///////////////////////////////////////////////
        		
        		String header = InputConvertorConstants.MS2PIP_HEADER_SPECID+"\t"+
        						InputConvertorConstants.MS2PIP_HEADER_MODIFICATIONS+"\t"+
        						InputConvertorConstants.MS2PIP_HEADER_PEPTIDE+"\t"+
        						InputConvertorConstants.MS2PIP_HEADER_CHARGE;
        		BW.append(header);
        		BW.newLine();
        		
        		/////////////////////////////////////////////////////////////////////////
        		
        		
        		ArrayList<MS2PIPRecord> records = new ArrayList<MS2PIPRecord>();
        		int discardedPRSMs = 0;
        		while((line = BR.readLine()) != null) {
        			String[] fields = line.split("\t");
        			
        			Peptide icPeptide	=	new Peptide(fields[icPeptideIdx], InputConvertorConstants.IC_CONVERTOR);
        			String infPeptide	=	fields[infPeptideIdx];
        			String charge		=	fields[chargeIdx];
        			
        			MS2PIPRecord record = new MS2PIPRecord();
        			
        			// IL determination from infPeptide
        			icPeptide.icPeptideToILDetermination(infPeptide);
        			// get modifications
        			// 1|Oxidation... something like that.
        			String modifications = MS2PIPRecord.getModifications(icPeptide);
        			
        			/****************************************************************************
        			 * #TODO #WARNING variable modifications
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
        			
        			record.wildPeptide = icPeptide.stripPeptide;
        			record.modifications = modifications;
        			record.charge = charge;
        			record.fullRecord = line;
        			
        			// add record
        			records.add(record);
        			record.idx = records.size();
        		}
        		
        		System.out.println("A total of PRSMs: "+records.size());
        		if(discardedPRSMs != 0) {
        			System.out.println("Note! "+discardedPRSMs+" PRSMs were discarded due to having unsupported modifications");
        		}
        		
        		// remove duplicated peptides + charge
        		Hashtable<String, MS2PIPRecord> hasPeptide = new Hashtable<String, MS2PIPRecord>();
        		for(int i=0; i<records.size(); i++) {
        			MS2PIPRecord record = records.get(i);
        			record.key = MS2PIPRecord.getMS2PIPKey(record.wildPeptide, record.modifications, record.charge);
        			MS2PIPRecord hasRecord = hasPeptide.get(record.key);
        			
        			if(hasRecord == null) {
        				hasPeptide.put(record.key, record);
        			}
        			
        		}
        		//// clear the records and insult again from hasPeptide
        		records.clear();
        		hasPeptide.forEach((p, record)->{
        			records.add(record);
        		});
        		// sort by index to make sure the order.
        		Collections.sort(records);
        		System.out.println("A total of peptides to be considered in MS2PIP: "+records.size());
        		
        		for(int i=0; i<records.size(); i++) {
        			MS2PIPRecord record = records.get(i);
        			BW.append(record.key).append("\t")
        			.append(record.modifications).append("\t")
        			.append(record.wildPeptide).append("\t")
        			.append(record.charge);
        			BW.newLine();
        		}
        		
        		BR.close();
        		BW.close();
        		
        		if(isDebugMode) {
        			BW = new BufferedWriter(new FileWriter(outputPath.replace(".input", ".debug")));
        			BW.append(header).append("\t").append(pXgHeader);
        			BW.newLine();
        			
        			for(int i=0; i<records.size(); i++) {
        				MS2PIPRecord record = records.get(i);
            			BW.append(record.key).append("\t")
            			.append(record.modifications).append("\t")
            			.append(record.wildPeptide).append("\t")
            			.append(record.charge).append("\t")
            			.append(record.fullRecord);
            			BW.newLine();
            		}
        			
        			BW.close();
        		}
        	}
        }
	}
	
	public void toAutoRTInputFormat (CommandLine cmd) throws IOException, ParseException {
		
		boolean isDebugMode = true;
		
		// Options //////////
		
		/**
		 * -i -P -I -R -S -p
		 * 
		 * -i pXg/
		 * -p Set01
		 * -S 3
		 * -R 10
		 * -P 12
		 * -I 20
		 * 
		 * For example>
		 * 
		 * target files:
		 * pXg/abc_Set01.pXg
		 * pXg/dfg_Set01.pXg
		 * pXg/ggg_set01.pXg // case sensitive
		 * ...
		 * 
		 * result files
		 * pXg/abc_Set01.autort.input
		 * pXg/dfg_Set01.autort.input
		 */
		/////////////////////
		
        String inputFile	= cmd.getOptionValue("i");
        String filePattern	= cmd.getOptionValue("p");
        int scoreIdx		= Integer.parseInt(cmd.getOptionValue("S"));
        int icRTIdx			= -1;
        int icPeptideIdx 	= -1;
        int infPeptideIdx	= -1;
        
        File[] files = new File(inputFile).listFiles();
        if(files == null) {
        	System.out.println("-i option should be a path of folder including .pXg files");
        	System.exit(1);
        }
        
        for(File file : files) {
        	if(file.getName().startsWith(".")) continue;
        	if(file.getName().endsWith(".pXg") && file.getName().contains(filePattern)) {
        		
        		BufferedReader BR = new BufferedReader(new FileReader(file));
        		System.out.println("Read "+file.getName());
        		String line = null;
        		String pXgHeader = BR.readLine(); // skip header
        		
        		// get index
        		String[] pXgHeaderSplit = pXgHeader.split("\t");
        		icRTIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_RT_FIELD_NAME);
        		icPeptideIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_PEPTIDE_FIELD_NAME);
        		infPeptideIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_INFERREND_PEPTIDE_FIELD_NAME);
        		
        		String outputPath = file.getAbsolutePath().replace(".pXg", ".autort.input");
        		BufferedWriter BW = new BufferedWriter(new FileWriter(outputPath));
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
        		
        		if(isDebugMode) {
        			BW = new BufferedWriter(new FileWriter(outputPath.replace(".input", ".debug")));
        			BW.append(header).append("\t").append(pXgHeader);
        			BW.newLine();
        			
        			for(int i=0; i<records.size(); i++) {
            			AutoRTRecord record = records.get(i);
            			BW.append(record.modifiedPeptide+"\t"+record.rt+"\t").append(record.fullRecord);
            			BW.newLine();
            		}
        			
        			BW.close();
        		}
        	}
        }
	}



}
