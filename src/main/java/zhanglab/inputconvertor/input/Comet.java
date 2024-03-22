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
import zhanglab.inputconvertor.data.ModificationTable;
import zhanglab.inputconvertor.data.Peptide;
import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.module.ToAutoRTInput;
import zhanglab.inputconvertor.module.ToMS2PIPInput;

public class Comet implements ToMS2PIPInput{
	public Comet () {}
	
	public void toMS2PIPInputFormat (CommandLine cmd) throws IOException, ParseException {
		boolean isDebugMode = true;
		
		// Options //////////
		
		/**
		 * -i -p -P -C
		 * 
		 * -i comet/
		 * -p Set01
		 * -P 12
		 * -C 2
		 * 
		 * For example>
		 * 
		 * target files:
		 * comet/abc_Set01.txt
		 * comet/dfg_Set01.txt
		 * comet/dfg_Set02.txt
		 * ...
		 * 
		 * result files
		 * comet/abc_Set01.ms2pip.input
		 * comet/dfg_Set01.ms2pip.input
		 */
		/////////////////////
		
        String inputFile	= cmd.getOptionValue("i");
        String filePattern	= cmd.getOptionValue("p");
        int chargeIdx		= Integer.parseInt(cmd.getOptionValue("C"));
        int modPeptideIdx 	= -1;
        
        File[] files = new File(inputFile).listFiles();
        if(files == null) {
        	System.out.println("-i option should be a path of folder including .txt files");
        	System.exit(1);
        }
        
        for(File file : files) {
        	if(file.getName().startsWith(".")) continue;
        	if(file.getName().endsWith(".txt") && file.getName().contains(filePattern)) {
        		
        		BufferedReader BR = new BufferedReader(new FileReader(file));
        		System.out.println("Read "+file.getName());
        		String line = null;
        		String cometHeader = BR.readLine(); // skip header
        		
        		String[] cometHeaderSplit = cometHeader.split("\t");
        		modPeptideIdx = InputConvertorConstants.getFieldIndex(cometHeaderSplit, InputConvertorConstants.COMET_PEPTIDE_FIELD_NAME);
        		
        		
        		String outputPath = file.getAbsolutePath().replace(".txt", ".ms2pip.input");
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
        		int discardedPSMs = 0;
        		while((line = BR.readLine()) != null) {
        			String[] fields = line.split("\t");
        			
        			Peptide peptide	=	new Peptide(fields[modPeptideIdx], InputConvertorConstants.COMET);
        			String charge		=	fields[chargeIdx];
        			
        			MS2PIPRecord record = new MS2PIPRecord();
        			
        			// get modifications
        			// 1|Oxidation... something like that.
        			String modifications = MS2PIPRecord.getModifications(peptide);
        			
        			/****************************************************************************
        			 * #TODO #WARNING variable modifications
        			 * Current version ignores other PTMs except for M_OXI.
        			 * Just drop the records!
        			 ****************************************************************************/
        			//check if there are other PTMs...
        			String checkPeptide = peptide.stripPeptide.replaceAll("[+-0123456789.*]", "");
        			if(!checkPeptide.equalsIgnoreCase(peptide.stripPeptide)) {
        				// drop the peptide
        				discardedPSMs++;
        				continue;
        			}
        			
        			record.wildPeptide = peptide.stripPeptide;
        			record.modifications = modifications;
        			record.charge = charge;
        			record.fullRecord = line;
        			
        			// add record
        			records.add(record);
        			record.idx = records.size();
        		}
        		
        		System.out.println("A total of PSMs: "+records.size());
        		if(discardedPSMs != 0) {
        			System.out.println("Note! "+discardedPSMs+" PSMs were discarded due to having unsupported modifications");
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
        			BW.append(header).append("\t").append(cometHeader);
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
		 * -i -P -I -R -S
		 * 
		 * -i comet/
		 * -p Set01
		 * -S 6
		 * -R 18
		 * -P 12
		 * 
		 * For example>
		 * 
		 * target files:
		 * comet/abc_Set01.txt
		 * comet/dfg_Set01.txt
		 * comet/dfg_Set02.txt
		 * ...
		 * 
		 * result files
		 * comet/abc_Set01.autort.input
		 * comet/dfg_Set01.autort.input
		 */
		/////////////////////
		
        String inputFile	= cmd.getOptionValue("i");
        String filePattern	= cmd.getOptionValue("p");
        int scoreIdx		= Integer.parseInt(cmd.getOptionValue("S"));
        int icRTIdx			= -1;
        int modPeptideIdx 	= -1;
        
        File[] files = new File(inputFile).listFiles();
        if(files == null) {
        	System.out.println("-i option should be a path of folder including .txt files");
        	System.exit(1);
        }
        
        for(File file : files) {
        	if(file.getName().startsWith(".")) continue;
        	if(file.getName().endsWith(".txt") && file.getName().contains(filePattern)) {
        		
        		BufferedReader BR = new BufferedReader(new FileReader(file));
        		System.out.println("Read "+file.getName());
        		String line = null;
        		String cometHeader = BR.readLine(); // skip header
        		
        		// index
        		String[] cometHeaderSplit = cometHeader.split("\t");
        		icRTIdx = InputConvertorConstants.getFieldIndex(cometHeaderSplit, InputConvertorConstants.COMET_RT_FIELD_NAME);
        		modPeptideIdx = InputConvertorConstants.getFieldIndex(cometHeaderSplit, InputConvertorConstants.COMET_PEPTIDE_FIELD_NAME);
        		
        		
        		String outputPath = file.getAbsolutePath().replace(".txt", ".autort.input");
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
        			Peptide peptide	=	new Peptide(fields[modPeptideIdx], InputConvertorConstants.COMET);
        			
        			AutoRTRecord record = new AutoRTRecord();
        			
        			record.score = score;
        			// sec to min
        			record.rt = (rt/60)+"";
        			
        			/****************************************************************************
        			 * TODO: variable modifications
        			 * Current version ignores other PTMs except for M_OXI.
        			 * Just drop the records!
        			 ****************************************************************************/
        			//check if there are other PTMs...
        			String checkPeptide = peptide.stripPeptide.replaceAll("[+-0123456789.*]", "");
        			if(!checkPeptide.equalsIgnoreCase(peptide.stripPeptide)) {
        				// drop the peptide
        				discardedPRSMs++;
        				continue;
        			}
        			
        			// convert IC_M_OXI to AUTORT_M_OXI
        			// e.g. ACM+15.995GG => AC1GG
        			peptide.toAutoRTModPeptide();
        			
        			record.modifiedPeptide = peptide.modPeptide;
        			record.fullRecord = line;
        			
        			// add record
        			records.add(record);
        			record.idx = records.size();
        		}
        		
        		System.out.println("A total of PSMs: "+records.size());
        		if(discardedPRSMs != 0) {
        			System.out.println("Note! "+discardedPRSMs+" PSMs were discarded due to having unsupported modifications");
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
        		System.out.println("Remain top-scored PSM per peptide ... ");
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
        			BW.append(header).append("\t").append(cometHeader);
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
