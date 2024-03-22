package zhanglab.inputconvertor.module;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;

import zhanglab.inputconvertor.env.InputConvertorConstants;

public class CompactDatabase {

	public static ArrayList<String> inputPaths = new ArrayList<String>();
	public static String outputFilePath;
	public static String decoyPrefix = "XXX";

	private CompactDatabase () {}
	
	public static void merge(String[] args) throws IOException {
		parseOptions(args);
		
		Hashtable<String, Boolean> checkDup = new Hashtable<String, Boolean>();
		
		File outputFile = new File(outputFilePath);
		
		BufferedWriter BW = new BufferedWriter(new FileWriter(outputFile));
		
		long totalSequence = 0;
		for(String fileName : inputPaths) {
			File file = new File(fileName);
			if(file.getName().startsWith(".")) continue;
			
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String line = null;
			
			String mark = "";
			while((line = BR.readLine()) != null) {
				if(line.startsWith(">")) {
					mark = line.substring(1).split("\\_")[0];
				} else {
					totalSequence++;
					String sequence = line;
					boolean isReference = mark.startsWith(InputConvertorConstants.REF_HEADER_ID) ? true : false;
					if(checkDup.get(sequence) == null || checkDup.get(sequence) == false) {
						checkDup.put(sequence, isReference);
					}
				}
			}
			BR.close();
		}
		
		Iterator<String> sequences = (Iterator<String>) checkDup.keys();
		// target
		int idx = 1;
		while(sequences.hasNext()) {
			String sequence = sequences.next();
			String mark = InputConvertorConstants.REF_HEADER_ID;
			
			if(checkDup.get(sequence) == false) {
				mark = InputConvertorConstants.NON_REF_HEADER_ID;
			}
			
			BW.append(">"+mark+"_"+idx++);
			BW.newLine();
			BW.append(sequence);
			BW.newLine();
		}
		
		// decoy
		sequences = (Iterator<String>) checkDup.keys();
		idx = 1;
		StringBuilder sb = new StringBuilder();
		while(sequences.hasNext()) {
			String sequence = sequences.next();
			String mark = InputConvertorConstants.REF_HEADER_ID;
			
			if(checkDup.get(sequence) == false) {
				mark = InputConvertorConstants.NON_REF_HEADER_ID;
			}
			
			BW.append(">"+decoyPrefix+"_"+mark+"_"+idx++);
			BW.newLine();
			BW.append(sb.append(sequence).reverse().toString());
			BW.newLine();
			sb.setLength(0);
		}
		
		System.out.println("Compaction DB: "+totalSequence +" >> "+checkDup.size());
		
		BW.close();
	}

	public static void parseOptions (String[] args) {
	    for(int i=0; i<args.length; i++) {
	    	if(args[i].equalsIgnoreCase("-o") || args[i].equalsIgnoreCase("--output")) {
	    		outputFilePath = args[++i];
	    	} else if(args[i].equalsIgnoreCase("-d") || args[i].equalsIgnoreCase("--decoy")) {
	    		decoyPrefix = args[++i];
	    	} else if(args[i].startsWith("-")){
	    		i++;
	    	} else {
	    		inputPaths.add(args[i]);
	    	}
	    }
	    
		System.out.println();
	}
}
