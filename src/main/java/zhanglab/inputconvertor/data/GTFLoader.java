package zhanglab.inputconvertor.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

public class GTFLoader {
	// gene id with version to transcripts
	public Hashtable<String, ArrayList<Transcript>> geneToTranscripts = new Hashtable<String, ArrayList<Transcript>>();
	
	public GTFLoader (File file) {
		try {
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String line = null;
			
			while((line = BR.readLine()) != null) {
				if(line.startsWith("#")) {
					continue;
				}
				
				String[] fields = line.split("\t");
				String feature = fields[2];
				
				if(feature.equalsIgnoreCase("exon")) {
					String chr = fields[0];
					int start = Integer.parseInt(fields[3]);
					int end = Integer.parseInt(fields[4]);
					String strand = fields[6];
					String[] attrs = fields[8].split("\\;");
					String geneId = getTag(attrs, "gene_id");
					String enstId = getTag(attrs, "transcript_id");
					
					Transcript transcript = getTranscript(geneId, enstId);
					transcript.strand = strand;
					transcript.chr = chr;
					transcript.exons.add(new Exon(start, end));
					
				}
			}
			
			BR.close();
			
			geneToTranscripts.forEach((g, ts)->{
				for(Transcript t : ts) {
					Collections.sort(t.exons);
				}
			});
		}catch(IOException ioe) {
			
		}
	}
	
	private Transcript getTranscript (String geneId, String enstId) {
		Transcript thisT = null;
		ArrayList<Transcript> transcripts = this.geneToTranscripts.get(geneId);
		if(transcripts == null) {
			transcripts = new ArrayList<Transcript>();
			thisT = new Transcript();
			thisT.tID = enstId;

			transcripts.add(thisT);
			this.geneToTranscripts.put(geneId, transcripts);
		} else {
			for(int i=0; i<transcripts.size(); i++) {
				Transcript thatT = transcripts.get(i);
				if(thatT.tID.equalsIgnoreCase(enstId)) {
					thisT = thatT;
					break;
				}
			}
			
			if(thisT == null) {
				thisT = new Transcript();
				thisT.tID = enstId;
				transcripts.add(thisT);
			}
		}
		
		return thisT;
	}
	
	private String getTag (String[] attrs, String tag) {
		String val = null;
		
		for(int i=0; i<attrs.length; i++) {
			if(attrs[i].contains(tag)) {
				val = attrs[i].trim().split("\\s")[1].replace("\"", "");
			}
		}
		
		return val;
	}
}
