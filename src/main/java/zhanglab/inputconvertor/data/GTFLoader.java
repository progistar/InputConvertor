package zhanglab.inputconvertor.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

public class GTFLoader {
	// gene id with version to transcripts
	public Hashtable<String, ArrayList<Transcript>> geneToTranscripts = new Hashtable<String, ArrayList<Transcript>>();
	public Hashtable<String, String> geneToGeneName = new Hashtable<String, String>();
	
	public void clear() {
		geneToTranscripts.clear();
	}
	
	public GTFLoader (File file) {
		System.out.println("Load GTF: "+file.getName());
		try {
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String line = null;
			
			Transcript transcript = null;
			while((line = BR.readLine()) != null) {
				if(line.startsWith("#")) {
					continue;
				}
				
				String[] fields = line.split("\t");
				String feature = fields[2];
				
				// #FORMAT
				//// StringTie or Reference
				if(feature.equalsIgnoreCase("exon") || feature.equalsIgnoreCase("cds") || feature.equalsIgnoreCase("transcript")) {
					String chr = fields[0];
					int start = Integer.parseInt(fields[3]);
					int end = Integer.parseInt(fields[4]);
					String strand = fields[6];
					String[] attrs = fields[8].split("\\;");
					
					// for gencode GTF
					String geneId = getTag(attrs, "gene_id");
					String enstId = getTag(attrs, "transcript_id");
					String geneName = getTag(attrs, "gene_name");
					
					// for stringtie GTF
					String refGeneId = getTag(attrs, "ref_gene_id");
					String refEnstId = getTag(attrs, "reference_id"); // base stringtie output
					if(refEnstId == null) {
						refEnstId = getTag(attrs, "cmp_ref"); // gff_comapre output
					}
					String refGeneName = getTag(attrs, "ref_gene_name");
					
					if(refGeneId != null) {
						geneId = refGeneId;
					}
					if(refEnstId != null) {
						enstId = refEnstId;
					}
					if(refGeneName != null) {
						geneName = refGeneName;
					}
					
					if(feature.equalsIgnoreCase("transcript")) {
						// enroll mapping table
						if(geneId != null && geneName != null) {
							geneToGeneName.put(geneId, geneName);
						}
						
						transcript = getTranscript(geneId, enstId);
						
						String fpkm = getTag(attrs, "FPKM");
						double expValue = Double.MAX_VALUE;
						if(fpkm != null) {
							expValue = Double.parseDouble(fpkm);
						}
						transcript.strand = strand;
						transcript.chr = chr;
						transcript.start = start+"";
						transcript.end = end+"";
						transcript.attrs = fields[8];
						transcript.FPKM = expValue;
						
						String classCode = getTag(attrs, "class_code");
						if(classCode == null) {
							transcript.classCode = "na";
						} else {
							transcript.classCode = classCode;
						}
					} else {
						if(feature.equalsIgnoreCase("exon")) {
							transcript.exons.add(new Exon(chr, start, end));
						} else if(feature.equalsIgnoreCase("cds")) {
							// CDS included
							transcript.isProteinCoding = true;
							transcript.cdss.add(new Exon(chr, start, end));
						}
					}
				}
				// #FORMAT
				//// CIRIquant
				else if(feature.equalsIgnoreCase("circRNA")) {
					String chr = fields[0];
					int start = Integer.parseInt(fields[3]);
					int end = Integer.parseInt(fields[4]);
					String strand = fields[6];
					String[] attrs = fields[8].split("\\;");
					String geneId = getTag(attrs, "gene_id");
					String enstId = getTag(attrs, "circ_id").replace("|", ",");
					
					// CIRIquant can have mutiple gene ids per record
					ArrayList<String> geneIds = new ArrayList<String>();
					if(geneId == null) {
						geneIds.add("intergenic");
					} else {
						String[] ids = geneId.split("\\,");
						for(String id : ids) {
							geneIds.add(id);
						}
					}
					
					for(String id : geneIds) {
						transcript = getTranscript(id, enstId);
						transcript.strand = strand;
						transcript.chr = chr;
						transcript.start = start+"";
						transcript.end = end+"";
						transcript.attrs = fields[8];
						transcript.exons.add(new Exon(chr, start, end));
					}
				}
			}
			
			BR.close();
			
			geneToTranscripts.forEach((g, ts)->{
				for(Transcript t : ts) {
					// System.out.println(t.tID+"\t"+g+"\t"+t.chr);
					t.setExons();
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
				String thisTag = attrs[i].trim().split("\\s")[0];
				if(thisTag.equalsIgnoreCase(tag)) {
					val = attrs[i].trim().split("\\s")[1].replace("\"", "");
				}
			}
		}
		
		return val;
	}
}
