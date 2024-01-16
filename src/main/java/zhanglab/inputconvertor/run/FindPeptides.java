package zhanglab.inputconvertor.run;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class FindPeptides {

	public static void main(String[] args) throws IOException {
		int geneNameIdx = 42;
		int eventIdx = 44;
		int mutationIdx = 36;
		
		File[] files = new File("/Users/seunghyukchoi/Documents/1_Projects/2023_Neoflow2/1_LUAD/LUAD_casanovo+pXg").listFiles();
		
		
		String header = null;
		for(File file :files) {
			if(file.getName().startsWith(".")) continue;
			if(!file.getName().endsWith(".ba")) continue;
			
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
				if(fields[geneNameIdx].equalsIgnoreCase("MAVS")) {
					System.out.println(sample+"\t"+line);
				}
				/*
				if(fields[eventIdx].contains("AS")) {
					System.out.println(sample+"\t"+line);
				}*/
				/*
				if(!fields[mutationIdx].contains("-")) {
					System.out.println(sample+"\t"+line);
				}*/
			}
			
			BR.close();
		}
	}
}
