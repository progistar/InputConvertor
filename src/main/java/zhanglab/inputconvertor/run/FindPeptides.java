package zhanglab.inputconvertor.run;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class FindPeptides {

	public static void main(String[] args) throws IOException {
		int geneNameIdx = 41;
		int eventIdx = 43;
		int mutationIdx = 35;
		
		File[] files = new File("/Users/seunghyukchoi/Documents/1_Projects/2023_Spliceosome/1.Search/0.Bruker").listFiles();
		
		
		String header = null;
		for(File file :files) {
			if(file.getName().startsWith(".")) continue;
			if(!file.getName().endsWith(".pXg")) continue;
			
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String line = null;
			
			if(header == null) {
				header = BR.readLine();
				System.out.println("Sample\t"+header);
			} else {
				BR.readLine();
			}
			
			String sample = file.getName().split("\\.")[0];
			
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				/*
				if(fields[geneNameIdx].equalsIgnoreCase("NADK")) {
					System.out.println(sample+"\t"+line);
				}*/
				
				if(fields[eventIdx].contains("AS")) {
					System.out.println(sample+"\t"+line);
				}
				/*
				if(!fields[mutationIdx].contains("-")) {
					System.out.println(sample+"\t"+line);
				}*/
			}
			
			BR.close();
		}
	}
}
