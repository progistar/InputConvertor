package zhanglab.inputconvertor.run;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;

import zhanglab.inputconvertor.data.FastaEntry;
import zhanglab.inputconvertor.data.FastaLoader;

public class TempRunSubsetByRPHM {

	public static void main(String[] args) throws IOException {
		
		String[] samples = {"C3L-06706", "C3L-06766", "C3L-06774", "C3L-06922", "C3L-07146"};
		
		File[] files = new File("/Volumes/Seunghyuk/CPTAC_STAD/CPTAC_STAD_Unmapped_DB").listFiles();
		
		for(File file : files) {
			
			if(file.getName().startsWith(".")) continue;
			if(!file.getName().endsWith(".fasta")) continue;
			
			boolean isIncluded = false;
			
			String sample = file.getName();
			
			for(int i=0; i<samples.length; i++) {
				if(sample.contains(samples[i])) {
					isIncluded = true;
				}
			}
			
			if(!isIncluded) {
				continue;
			}
			
			System.out.println(sample);
			

			FastaLoader fasta = new FastaLoader(file);
			
			for(FastaEntry fastaEntry : fasta.entries) {
				fastaEntry.geneId = fastaEntry.originHeader.split("\\|")[2].split("\\s")[0];
			}
			
			// sort by geneId (=RPHM)
			Collections.sort(fasta.entries, new Comparator<FastaEntry>() {
				@Override
				public int compare(FastaEntry o1, FastaEntry o2) {
					double rphm1 = Double.parseDouble(o1.geneId);
					double rphm2 = Double.parseDouble(o2.geneId);
							
					if(rphm1 > rphm2) {
						return -1;
					} else if(rphm1 < rphm2) {
						return 1;
					}
							
					return 0;
				}
			});
			
			
			BufferedWriter BW = new BufferedWriter(new FileWriter(file.getAbsolutePath().replace(".fasta", ".top100k.fasta")));
		        
			int cnt = 0;
	        for(FastaEntry entry : fasta.entries) {
	        	cnt ++;
	        	if(cnt > 100000) {
	        		break;
	        	}
	        	BW.append(">"+entry.originHeader);
	        	BW.newLine();
	        	BW.append(entry.sequence);
	        	BW.newLine();
	        }
	        
	        BW.close();
		}
		
		
	}
	
}
