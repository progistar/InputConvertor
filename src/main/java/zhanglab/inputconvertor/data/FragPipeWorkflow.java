package zhanglab.inputconvertor.data;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class FragPipeWorkflow {

	private ArrayList<String> contents = new ArrayList<String>();
	
	public FragPipeWorkflow (File file) throws IOException {
		System.out.println("Read "+file.getName());
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		while((line = BR.readLine()) != null) {
			contents.add(line);
		}
		
		BR.close();
	}
	
	public void setParam (String parameterName, String value) {
		for(int i=0; i<contents.size(); i++) {
			String content = contents.get(i);
			
			if(content.contains("=")) {
				String param = content.split("\\=")[0];
				if(param.equalsIgnoreCase(parameterName)) {
					param += "=" + value;
					contents.set(i, param);
				}
			}
		}
	}
	
	public void write (File outputFile) throws IOException {
		System.out.println("Write "+outputFile.getName());
		BufferedWriter BW = new BufferedWriter(new FileWriter(outputFile));
		
		for(String line : contents) {
			BW.append(line);
			BW.newLine();
		}
		
		BW.close();
	}
}
