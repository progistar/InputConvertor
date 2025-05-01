package zhanglab.inputconvertor.run;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class TitleReformat {

	public static boolean ADD_SCANS = false;
	public static boolean REMOVE_UNKNOWN_CHARGE = true;
	
	public static void main(String[] args) throws IOException {
		File[] files = new File("/Volumes/Papers/2023_Spliceosome/0.Zips/SUM159_Bruker_MGF").listFiles();
		
		for(File file : files) {
			if(file.getName().endsWith(".mgf") && !file.getName().startsWith(".")) {
				System.out.println(file.getName());
				BufferedReader BR = new BufferedReader(new FileReader(file));
				BufferedWriter BW = new BufferedWriter(new FileWriter(file.getAbsolutePath().replace(".mgf", ".rmUnknownCharge.mgf")));
				BufferedWriter BW_unkonwn = new BufferedWriter(new FileWriter(file.getAbsolutePath().replace(".mgf", ".unknownCharge.mgf")));
				String line = null;
				double known = 0;
				double unknown = 0;
				
				ArrayList<String> spectrum = new ArrayList<String>();
				boolean isCharged = false;
				
				while((line = BR.readLine()) != null) {
					
					if(line.startsWith("BEGIN")) {
						spectrum = new ArrayList<String>();
						isCharged = false;
					}
					
					if(line.startsWith("TITLE")) {
						line = line.split("\\s")[0];
						String scanNum = line.split("\\.")[2];
						spectrum.add(line);
						if(ADD_SCANS) {
							spectrum.add("SCANS="+scanNum);	
						}
					} else if (line.startsWith("CHARGE=") ) {
						isCharged = true;
						spectrum.add(line);
					} else {
						spectrum.add(line);
					}
					
					
					if(line.startsWith("END")) {
						if(REMOVE_UNKNOWN_CHARGE && !isCharged) {
							unknown++;
							for(String str : spectrum) {
								BW_unkonwn.append(str);
								BW_unkonwn.newLine();
							}
						} else {
							known++;
							for(String str : spectrum) {
								BW.append(str);
								BW.newLine();
							}
						}
					}
				}
				BW_unkonwn.close();
				BW.close();
				BR.close();
				
				System.out.println(known/(unknown+known));
			}
		}
	}
}
