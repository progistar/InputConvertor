package zhanglab.inputconvertor.run;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import zhanglab.inputconvertor.env.InputConvertorConstants;

public class CheckFragmented {

	public static void main(String[] args) throws IOException {
		String filePath = "/Users/seunghyukchoi/Documents/1_Projects/2023_Spliceosome/1.Search/0.JMChoi_MSFragger/OT_Astral_Fragpipe_metadata/Fragpipe_metadata/8800_1/peptide.tsv";
		File[] files = new File(filePath).listFiles();
		if(files == null) {
			files = new File[1];
			files[0] = new File(filePath);
		}
		
		String targetFilePrefix = "Bruker_DMSO_Rep";
//		String targetFilePrefix = "BCM_DMSO_Rep";
//		String targetFilePrefix = "BCM_SF3B1_dTAG3_Rep";
		Hashtable<String, String>  targetSevenMers = new Hashtable<String, String>();
		Hashtable<String, String>  targetOtherMers = new Hashtable<String, String>();
		for(File file :files) {
			if(file.getName().startsWith(".")) continue;
			//if(!file.getName().endsWith(".fdr")) continue;
			
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String line = null;
			
			String[] header = BR.readLine().split("\t");
			int infPeptIdx = InputConvertorConstants.getFieldIndex(header, "Peptide");
			
			Hashtable<String, String>  sevenMers = new Hashtable<String, String>();
			Hashtable<String, String>  otherMers = new Hashtable<String, String>();
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String peptide = fields[infPeptIdx];
				if(peptide.length() == 7) {
					sevenMers.put(peptide, "non-overlap");
					if(file.getName().startsWith(targetFilePrefix)) {
						targetSevenMers.put(peptide, "non-overlap");
					}
				} else if (peptide.length() >= 8 && peptide.length() <= 15){
					otherMers.put(peptide, "");
					if(file.getName().startsWith(targetFilePrefix)) {
						targetOtherMers.put(peptide, "");
					}
				}
			}
			BR.close();
			
			
			printOverlapStat(sevenMers, otherMers, file.getName());
		}
		printOverlapStat(targetSevenMers, targetOtherMers, targetFilePrefix);
	}
	
	public static void printOverlapStat (Hashtable<String, String> sevenMers, Hashtable<String, String> otherMers, String name) {
		ArrayList<String> sevens = new ArrayList<String>();
		ArrayList<String> others = new ArrayList<String>();
		sevenMers.forEach((p, t)->{
			sevens.add(p);
		});
		otherMers.forEach((p, t)->{
			others.add(p);
		});
		// do
		for(String p : sevens) {
			for(String t : others) {
				if(t.contains(p)) {
					sevenMers.put(p, "overlap");
				}
			}
		}
		
		// 0 is non-overlap
		// 1 is overlap
		double[] overlapCount = new double[2];
		sevenMers.forEach((p, t)->{
			if(t.equalsIgnoreCase("overlap")) {
				overlapCount[1]++;
			} else {
				overlapCount[0]++;
			}
		});
		double overlapRatio = overlapCount[1] / (overlapCount[0]+overlapCount[1]);
		System.out.println(name+"\t"+overlapCount[0]+"\t"+overlapCount[1]+"\t"+overlapRatio);
	}
}
