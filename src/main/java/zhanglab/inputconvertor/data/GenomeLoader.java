package zhanglab.inputconvertor.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

public class GenomeLoader {

	private Hashtable<String, StringBuilder> genomeMap = new Hashtable<String, StringBuilder>();
	
	public GenomeLoader (File file) {
		// load genome
		try {
			
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String line = null;
			
			StringBuilder sequence = null;
			while((line = BR.readLine()) != null) {
				if(line.startsWith(">")) {
					String chr = line.split("\\s")[0].substring(1);
					System.out.println("Read "+chr);
					sequence = new StringBuilder();
					genomeMap.put(chr, sequence);
				} else {
					sequence.append(line.toUpperCase());
				}
			}
			
			BR.close();
			
		}catch(IOException ioe) {
			
		}
	}
	
	/**
	 * zero-based, [start, end)
	 * 
	 * @param chr
	 * @param start
	 * @param end
	 * @return
	 */
	public String getSequence(String chr, int start, int end) {
		StringBuilder sequence = genomeMap.get(chr);
		if(sequence == null) {
			return null;
		}
		
		return sequence.substring(start, end);
	}
}
