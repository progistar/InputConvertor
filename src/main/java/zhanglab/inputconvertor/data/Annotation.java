package zhanglab.inputconvertor.data;

import java.util.ArrayList;
import java.util.Hashtable;

import org.ahocorasick.trie.Emit;

public class Annotation {

	
	public Hashtable<String, Hashtable<String, GenomicInformation>> map = new Hashtable<String, Hashtable<String, GenomicInformation>>();
	
	public void putAnnotation (Emit emit, String header, String fullSequence, boolean isContam) {
		Hashtable<String, GenomicInformation> genomicInformations = map.get(emit.getKeyword());
		if(genomicInformations == null) {
			genomicInformations = new Hashtable<String, GenomicInformation>();
			map.put(emit.getKeyword(), genomicInformations);
		}
		
		// ensure that
		// if a genomic location is overlapped to previous information, the better annotation is selected.
		GenomicInformation genomicInformation = new GenomicInformation(header, emit, fullSequence, isContam);
		GenomicInformation prevInformation = genomicInformations.get(genomicInformation.getKey());
		if(prevInformation == null) {
			genomicInformations.put(genomicInformation.getKey(), genomicInformation);
		} else {
			if(genomicInformation.penalty < prevInformation.penalty) {
				genomicInformations.put(genomicInformation.getKey(), genomicInformation);
			}
		}
	}
	
	public void calRepresentAnnotation () {
		map.forEach((peptide, gInfoMap)->{
			ArrayList<String> filterLocations = new ArrayList<String>();
			int[] minPenalty = new int[1];
			minPenalty[0] = Integer.MAX_VALUE;
			gInfoMap.forEach((location, gInfo)->{
				minPenalty[0] = Math.min(minPenalty[0], gInfo.penalty);
			});
			
			gInfoMap.forEach((location, gInfo)->{
				if(minPenalty[0] < gInfo.penalty) {
					filterLocations.add(location);
				}
			});
			
			for(String location : filterLocations) {
				gInfoMap.remove(location);
			}
		});
	}
	
	//"GenomicLoci\tGenomicLociCount\tGeneName\tGeneNameCount\tMutationCount\tCategory";
	public ArrayList<String> getModifiedRecordList (String fullRecord, String peptide) {
		ArrayList<String> modifiedList = new ArrayList<String>();
		Hashtable<String, GenomicInformation> gMap = map.get(peptide.replace("I", "L"));
		if(gMap == null) {
			System.out.println(peptide + " is out of the database");
			modifiedList.add("-\t0\t-\t0\tError");
		} else {
			StringBuilder meta = new StringBuilder();
			Hashtable<String, String> geneNames = new Hashtable<String, String>();
			// count gene name
			gMap.forEach((key, gInfo)->{
				String gName = gInfo.geneName;
				geneNames.put(gName, "");
			});
			
			final int genomicLociCount = gMap.size();
			final int geneNameCount = geneNames.size();
			
			gMap.forEach((key, gInfo)->{
				String genomicLoci = gInfo.location;
				String gName = gInfo.geneName;
				String category = gInfo.from;
				int mutCount = gInfo.mutCnt;
				
				meta.setLength(0);
				meta.append(fullRecord.replace(peptide, gInfo.ilPeptide)).append("\t").append(genomicLoci).append("\t")
				.append(genomicLociCount).append("\t")
				.append(gName).append("\t")
				.append(geneNameCount).append("\t")
				.append(mutCount).append("\t")
				.append(category);
				
				modifiedList.add(meta.toString());
			});
		}
		
		
		return modifiedList;
		
	}
}
