package zhanglab.inputconvertor.run;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

class PercentRecord implements Comparable<PercentRecord>{
	double score;
	double pScore;
	String record;
	@Override
	public int compareTo(PercentRecord o) {
		if(this.score < o.score) {
			return 1;
		} else if(this.score > o.score) {
			return -1;
		}
		return 0;
	}
}

public class RunPercentScore {


	public static void main(String[] args) throws IOException {
		File[] files = new File("/Users/gistar/Documents/ZhangLab/2023_Immunopeptidomics_LUAD/Test/4sample_test_pin").listFiles();
		
		for(File file : files) {
			if(file.getName().startsWith(".")) continue;
			if(!file.getName().endsWith(".pin")) continue;
			
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String outputPIN = file.getAbsolutePath().replace(".pin",".norm.pin");
			BufferedWriter BW = new BufferedWriter(new FileWriter(outputPIN));
			String line = null;
			
			BW.append(BR.readLine());
			BW.newLine();
			ArrayList<PercentRecord> records = new ArrayList<PercentRecord>();
			while((line = BR.readLine()) != null) {
				PercentRecord pRecord = new PercentRecord();
				pRecord.record = line;
				pRecord.score = Double.parseDouble(line.split("\t")[3]);
				records.add(pRecord);
			}
			
			Collections.sort(records);
			
			double size = records.size();
			for(int i=0; i<records.size(); i++) {
				PercentRecord pRecord = records.get(i);
				boolean isTie = false;
				if(i > 0) {
					PercentRecord pRecordPrev = records.get(i-1);
					if(pRecord.score == pRecordPrev.score) {
						isTie = true;
					}
				}
				if(isTie) {
					pRecord.pScore = records.get(i-1).pScore;
				} else {
					pRecord.pScore = (size-i) / size * 100;
				}
				
				
				String[] fields = pRecord.record.split("\t");
				fields[3] = pRecord.pScore+"";
				
				BW.append(fields[0]);
				for(int j=1; j<fields.length; j++) {
					BW.append("\t"+fields[j]);
				}
				BW.newLine();
			}
			
			BR.close();
			BW.close();
		}
	}
}
