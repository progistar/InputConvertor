package zhanglab.inputconvertor.run;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import zhanglab.inputconvertor.data.FastaEntry;
import zhanglab.inputconvertor.data.FastaLoader;
import zhanglab.inputconvertor.env.InputConvertorConstants;

public class MergeDatabase {

	public static void main(String[] args) throws IOException {
		File[] files = new File("/Users/seunghyukchoi/Documents/1_Projects/2023_Spliceosome/2.iRNA").listFiles();
		
		String[] sampleNames = {"8800", "DMSO", "FKBP-SF3B1_DMSO", "FKBP-SF3B1_dTag13"};
		String sufix = ".combined.fasta";
		for(String sampleName :sampleNames) {
			System.out.println("Target sample: "+sampleName);
			BufferedWriter BW = new BufferedWriter(new FileWriter(sampleName+sufix));
			for(File file : files) {
				if(file.getName().startsWith(".")) continue;
				if(file.getName().endsWith(".fasta") && file.getName().startsWith(sampleName)) {
					System.out.println(file.getName()+" is targeted");
					
					BufferedReader BR = new BufferedReader(new FileReader(file));
					
					String line = null;
					
					while((line = BR.readLine()) != null) {
						if(line.startsWith(">")) {
							// add file information
							line += "\t"+file.getName();
						}
						BW.append(line);
						BW.newLine();
					}
					
					BR.close();
				}
			}
			BW.close();
		}
		
		FastaLoader references = new FastaLoader(new File("/Users/seunghyukchoi/Documents/_resources/_databases/gencode.v42.pc_translations.fa"));
		FastaLoader contaminants = new FastaLoader(new File("/Users/seunghyukchoi/Documents/_resources/_databases/maxquant_contaminants.fasta"));
		
		for(String sampleName : sampleNames) {
			String fileName = sampleName+sufix;
			String compFileName = sampleName+".compressed"+sufix;
			Hashtable<String, Integer> isDuplicated = new Hashtable<String, Integer>();
			FastaLoader fl = new FastaLoader(new File(fileName));
			
			BufferedWriter BW = new BufferedWriter(new FileWriter(compFileName));
			int stIdx = 0;
			int irIdx = 0;
			int fgIdx = 0;
			int ciIdx = 0;
			ArrayList<FastaEntry> uniqueEntries = new ArrayList<FastaEntry>();
			Hashtable<String, String> dbFileCnt = new Hashtable<String, String>();
			for(FastaEntry f : fl.entries) {
				String dbFileName = f.originHeader.split("\t")[1];
				dbFileCnt.put(dbFileName, "");
				
				if(isDuplicated.get(f.sequence) == null) {
					isDuplicated.put(f.sequence, 1);
					uniqueEntries.add(f);
					String[] fields = f.originHeader.split("\\|");
					
					if(f.originHeader.startsWith(InputConvertorConstants.ARRIBA_HEADER_ID)) {
						fgIdx++;
						fields[0] = InputConvertorConstants.ARRIBA_HEADER_ID+fgIdx;
					} else if(f.originHeader.startsWith(InputConvertorConstants.STRINGTIE_HEADER_ID)) {
						stIdx++;
						fields[0] = InputConvertorConstants.STRINGTIE_HEADER_ID+stIdx;
					} else if(f.originHeader.startsWith(InputConvertorConstants.IRFINDER_HEADER_ID)) {
						irIdx++;
						fields[0] = InputConvertorConstants.IRFINDER_HEADER_ID+irIdx;
					} else if(f.originHeader.startsWith(InputConvertorConstants.CIRIQUANT_HEADER_ID)) {
						ciIdx++;
						fields[0] = InputConvertorConstants.CIRIQUANT_HEADER_ID+ciIdx;
					}
					f.originHeader = fields[0];
				} else {
					Integer cnt = isDuplicated.get(f.sequence);
					isDuplicated.put(f.sequence, cnt+1);
				}
			}
			
			int totalFiles = dbFileCnt.size();
			int idx = 1;
			for(FastaEntry f : references.entries) {
				BW.append(">").append(InputConvertorConstants.REF_HEADER_ID+idx);
				BW.newLine();
				BW.append(f.sequence);
				BW.newLine();
				idx++;
			}
			
			 idx = 1;
			for(FastaEntry f : contaminants.entries) {
				BW.append(">").append("CT"+idx);
				BW.newLine();
				BW.append(f.sequence);
				BW.newLine();
				idx++;
			}
			
			for(FastaEntry f : uniqueEntries) {
				Integer cnt = isDuplicated.get(f.sequence);
				BW.append(">").append(f.originHeader).append("_").append(cnt+"/"+totalFiles);
				BW.newLine();
				BW.append(f.sequence);
				BW.newLine();
			}
			
			BW.close();
		}
	}
} 
