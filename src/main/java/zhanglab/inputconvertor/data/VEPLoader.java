package zhanglab.inputconvertor.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.SortedMap;
import java.util.TreeMap;

import zhanglab.inputconvertor.env.InputConvertorConstants;

public class VEPLoader {
	
	public Hashtable<String, TreeMap<Integer, ArrayList<Mutation>>> positionalMap = new Hashtable<String, TreeMap<Integer, ArrayList<Mutation>>>();
	public VEPLoader (File file, boolean isSomatic) throws IOException {
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		String[] header = BR.readLine().split("\t");
		
		int chrIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.VEP_CHR_FIELD_NAME);
		int posIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.VEP_POS_FIELD_NAME);
		int refIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.VEP_REF_FIELD_NAME);
		int altIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.VEP_ALT_FIELD_NAME);
		
		int[] totalOfMutations = new int[3];
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String chr = fields[chrIdx];
			int pos = Integer.parseInt(fields[posIdx]);
			String ref = fields[refIdx];
			String alt = fields[altIdx];
			
			//**
			// VEP represents an alteration one by one.
			
			// mutation
			Mutation mutation = new Mutation();
			mutation.altSeq = alt; mutation.refSeq = ref; mutation.pos = pos; mutation.chr = chr;
			mutation.setMutationStatus();
			
			TreeMap<Integer, ArrayList<Mutation>> map = positionalMap.get(chr);
			if(map == null) {
				map = new TreeMap<Integer, ArrayList<Mutation>>();
				positionalMap.put(chr, map);
			}
			
			ArrayList<Mutation> mutations = map.get(mutation.pos);
			if(mutations == null) {
				mutations = new ArrayList<Mutation>();
				map.put(mutation.pos, mutations);
			}
			
			mutations.add(mutation);
			
			totalOfMutations[mutation.type]++;
		}
		
		BR.close();
		
		System.out.println("A total of mutations: "+ (totalOfMutations[InputConvertorConstants.SNP] + totalOfMutations[InputConvertorConstants.INS] + totalOfMutations[InputConvertorConstants.DEL]));
		System.out.println("SNPs: "+totalOfMutations[InputConvertorConstants.SNP]);
		System.out.println("INSs: "+totalOfMutations[InputConvertorConstants.INS]);
		System.out.println("DELs: "+totalOfMutations[InputConvertorConstants.DEL]);
	}
	
	/**
	 * zero-based <br>
	 * [start, end)
	 * @param start
	 * @param end
	 * @return
	 */
	public ArrayList<Mutation> getMutationByRange (String chr, int start, int end) {
		ArrayList<Mutation> mutations = new ArrayList<Mutation>();
		if(positionalMap.get(chr) != null) {
			SortedMap<Integer, ArrayList<Mutation>> sortedMap = positionalMap.get(chr).subMap(start, end);
			if(sortedMap != null) {
				sortedMap.forEach((pos, m)->{
					mutations.addAll(m);
				});
			}
		}
		
		
		return mutations;
	}
}
