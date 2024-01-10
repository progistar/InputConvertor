package zhanglab.inputconvertor.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import zhanglab.inputconvertor.env.InputConvertorConstants;

public class IRFinder {

	public IRFinder (File file, GTFLoader referenceGTF) throws IOException {
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String header = BR.readLine();
		String[] fields = header.split("\t");
		
		// #HEADER
		int chrIdx = InputConvertorConstants.getFieldIndex(fields, "Chr");
		int startIdx = InputConvertorConstants.getFieldIndex(fields, "Start");
		int endIdx = InputConvertorConstants.getFieldIndex(fields, "End");
		int strandIdx = InputConvertorConstants.getFieldIndex(fields, "Strand");
		int nameIdx = InputConvertorConstants.getFieldIndex(fields, "Name");
		int warningIdx = InputConvertorConstants.getFieldIndex(fields, "Warnings");
		//////////////////////////////////////// END OF #HEADER /////////////////
		
		String line = null;
		
		while((line = BR.readLine()) != null) {
			fields = line.split("\t");
			String chr = fields[chrIdx];
			String start = fields[startIdx];
			String end = fields[endIdx];
			String name = fields[nameIdx];
			String strand = fields[strandIdx];
			String warnings = fields[warningIdx];
			
			
			
		}
		
		
		BR.close();
	}
}
