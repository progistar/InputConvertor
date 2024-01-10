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
import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.data.Peptide;
import zhanglab.inputconvertor.env.InputConvertorConstants;

class TDRecord implements Comparable<TDRecord>{
	String record;
	double score;
	int label;
	double fdr;
	
	@Override
	public int compareTo(TDRecord o) {
		if(this.score < o.score) {
			return 1;
		} else if(this.score > o.score) {
			return -1;
		} else if(this.label > o.label) {
			return -1;
		} else if(this.label < o.label) {
			return 1;
		}
		return 0;
	}
	
	public String toString() {
		return record +"\t" + score;
	}
}

public class TargetDecoyAnalysis {

	public void doFDR (CommandLine cmd) throws IOException, ParseException {
		String inputFile = cmd.getOptionValue("i");
        String inputPattern = cmd.getOptionValue("p");
        double fdr = Double.parseDouble(cmd.getOptionValue("d"));
        
        File[] files = new File(inputFile).listFiles();
        if(files == null) {
        	System.out.println("-i option should be a path of folder including .pXg.feat files");
        	System.exit(1);
        }
        
        File targetFile = null;
        File decoyFile = null;
        File featFile = null;
        for(File file : files) {
        	if(file.getName().startsWith(".")) continue;
        	if(file.getName().endsWith(".pXg.feat") && file.getName().contains(inputPattern)) {
        		featFile = file;
        		System.out.println(featFile.getName());
        	} else if(file.getName().endsWith(".target.psm") && file.getName().contains(inputPattern)) {
        		targetFile = file;
        		System.out.println(targetFile.getName());
        	} else if(file.getName().endsWith(".decoy.psm") && file.getName().contains(inputPattern)) {
        		decoyFile = file;
        		System.out.println(decoyFile.getName());
        	}
        }
        
        Hashtable<String, String[]> targetPSMs = getHashTableFromPSM(targetFile);
        Hashtable<String, String[]> decoyPSMs = getHashTableFromPSM(decoyFile);
        
        
        BufferedReader BR = new BufferedReader(new FileReader(featFile));
        BufferedWriter BW = new BufferedWriter(new FileWriter(featFile.getAbsolutePath().replace(".feat", ".feat.fdr")));
        String line = null;
        String header = BR.readLine();
        
        // build header ////////////////////////////////////////////////
        BW.append(header).append("\t")
        .append("percoaltor_score");
        BW.newLine();
        ////////////////////////////////////// End of building header /
        
        String[] headerSplit = header.split("\t");
        int specIdIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.IC_SPEC_ID_FEILD_NAME);
        int labelIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.IC_LABEL_FEILD_NAME);
        int genomicIdIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.IC_GENOMIC_ID_FEILD_NAME);
        int isCanonicalIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.IC_IS_CANONICAL_FEILD_NAME);
        int icPeptideIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.IC_PEPTIDE_FIELD_NAME);
        int inferredPeptideIdx = InputConvertorConstants.getFieldIndex(headerSplit, InputConvertorConstants.IC_INFERREND_PEPTIDE_FIELD_NAME);
        
        
        ArrayList<TDRecord> cRecords = new ArrayList<TDRecord>();
        ArrayList<TDRecord> ncRecords = new ArrayList<TDRecord>();
        Hashtable<String, String> duplications = new Hashtable<String, String>();
        while((line = BR.readLine()) != null) {
        	String[] fields = line.split("\t");
        	String specId = fields[specIdIdx];
        	String genomicId = fields[genomicIdIdx];
        	String isCanonical = fields[isCanonicalIdx];
        	String label = fields[labelIdx];
        	Peptide peptide = new Peptide(fields[icPeptideIdx], InputConvertorConstants.IC_CONVERTOR);
        	peptide.icPeptideToILDetermination(fields[inferredPeptideIdx]);
        	
        	// key: specId + peptide.modPeptide
        	// this is because same genomic id and spec id can be occurred when there are modified peptides (site-localization problem).
        	// Casanovo sometimes distinguishes I/L... 
        	String key = specId+"_"+peptide.modPeptide;
        	
        	String[] psm = targetPSMs.get(key) != null ? targetPSMs.get(key) : decoyPSMs.get(key);
        	if(psm == null) {
        		continue;
        	}
        	
        	String genomicIdInPSM = psm[0].replace("XXX_", "");
        	
        	if(genomicIdInPSM.equalsIgnoreCase(genomicId)) {
        		TDRecord record = new TDRecord();
        		record.record = line;
        		record.score = Double.parseDouble(psm[1]);
        		record.label = Integer.parseInt(label);
        		
        		if(isCanonical.equalsIgnoreCase("true")) {
        			cRecords.add(record);
        		} else {
        			ncRecords.add(record);
        		}
        	}
        }
        ArrayList<TDRecord> passTargetList = new ArrayList<TDRecord>();
        passTargetList.addAll(getFDR(cRecords, fdr));
        passTargetList.addAll(getFDR(ncRecords, fdr));
        
        for(int i=0; i<passTargetList.size(); i++) {
        	TDRecord record = passTargetList.get(i);
        	
        	String specId = record.record.split("\t")[specIdIdx];
        	if(duplications.get(specId) == null) {
        		duplications.put(specId, specId);
        		BW.append(record.toString());
            	BW.newLine();
        	} else {
        		System.out.println("Duplicated records: "+specId);
        	}
        }
        
        BW.close();
        BR.close();
        
        
        
	}
	
	private ArrayList<TDRecord> getFDR (ArrayList<TDRecord> records, double fdr) {
		Collections.sort(records);
		
		double tCount = 0;
		double dCount = 0;
		int lastIdx   = 0;
		for(int i=0; i<records.size(); i++) {
			TDRecord record = records.get(i);
			
			if(record.label > 0) {
				tCount++;
			} else {
				dCount++;
			}
			
			if(tCount == 0) {
				record.fdr = 1;
			} else {
				record.fdr = dCount/tCount;
				
				if(record.fdr < fdr) {
					lastIdx = i;
				}
			}
		}
		
		ArrayList<TDRecord> passTargetList = new ArrayList<TDRecord>();
		for(int i=0; i<=lastIdx; i++) {
			TDRecord record = records.get(i);
			if(record.label > 0) {
				passTargetList.add(record);
			}
		}
		
		return passTargetList;
	}
	
	
	
	private Hashtable<String, String[]> getHashTableFromPSM (File file) throws IOException {
		Hashtable<String, String[]> table = new Hashtable<String, String[]>();
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		String[] header = BR.readLine().split("\t");
		int psmIdIdx = InputConvertorConstants.getFieldIndex(header, "PSMId");
		int scoreIdx = InputConvertorConstants.getFieldIndex(header, "score");
		int proteinIdIdx = InputConvertorConstants.getFieldIndex(header, "proteinIds");
		int peptideIdx = InputConvertorConstants.getFieldIndex(header, "peptide");
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String[] record = {fields[proteinIdIdx], fields[scoreIdx]};
			String key = fields[psmIdIdx]+"_"+fields[peptideIdx];
			table.put(key, record);
		}
		
		BR.close();
		
		return table;
	}
}
