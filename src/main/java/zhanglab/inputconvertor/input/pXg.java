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

import zhanglab.inputconvertor.data.DeepLCRecord;
import zhanglab.inputconvertor.data.MS2PIPRecord;
import zhanglab.inputconvertor.data.NetMHCpanResult;
import zhanglab.inputconvertor.data.Peptide;
import zhanglab.inputconvertor.data.Spectra;
import zhanglab.inputconvertor.data.Spectrum;
import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.function.CalculateSA;
import zhanglab.inputconvertor.function.NetMHCpanParser;
import zhanglab.inputconvertor.module.NetMHCpan;
import zhanglab.inputconvertor.module.ToAutoRTInput;
import zhanglab.inputconvertor.module.ToFeatures;
import zhanglab.inputconvertor.module.ToMS2PIPInput;


public class pXg {
	public pXg () {}
	
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
        		
        		infPeptideIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.IC_INFERRED_PEPTIDE_FIELD_NAME);
        		
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

	/**
	 * @deprecated
	 * 
	 * @param cmd
	 * @throws IOException
	 * @throws ParseException
	 */
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
        		
        		infPeptideIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.IC_INFERRED_PEPTIDE_FIELD_NAME);
        		
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
	
	/**
	 * @deprecated since we decided to use Prosit model.
	 * 
	 * */
	public void toMS2PIPInputFormat (CommandLine cmd) throws IOException, ParseException {
		boolean isDebugMode = false;
		
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
        int chargeIdx		= -1;
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
        		infPeptideIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_INFERRED_PEPTIDE_FIELD_NAME);
        		chargeIdx = InputConvertorConstants.getFieldIndex(pXgHeaderSplit, InputConvertorConstants.IC_CHARGE_FIELD_NAME);
        		
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
}