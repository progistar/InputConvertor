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
	public boolean isSomatic = false;
	
	public VEPLoader (File file, boolean isSomatic) throws IOException {
		this.isSomatic = isSomatic;
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		String[] header = BR.readLine().split("\t");
		
		int locIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.VEP_LOC_FIELD_NAME);
		int refIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.VEP_REF_FIELD_NAME);
		int altIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.VEP_ALT_FIELD_NAME);
		
		int[] totalOfMutations = new int[3];
		Hashtable<String, String> removeDuplication = new Hashtable<String, String>();
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String location = fields[locIdx];
			String chr = location.split("\\:")[0];
			if(!chr.startsWith("chr")) {
				chr = "chr"+chr;
			}
			String pos = location.split("\\:")[1];
			int startPos = Integer.parseInt(pos.split("\\-")[0]);
			int endPos = Integer.parseInt(pos.split("\\-")[1]);
			String refs = fields[refIdx].replace("-", "");
			String alts = fields[altIdx].replace("-", "");
			int refLen = refs.length();
			int altLen = alts.length();
			
			String key = chr+":"+pos+refs+">"+alts;
			if(removeDuplication.get(key) != null) {
				continue;
			}
			removeDuplication.put(key, "");
			
			//**
			// VEP can contain di/tri something like that
			
			// SNP / DNP / TNP ...
			if(refLen == altLen) {
				int idx = 0;
				for(int i=startPos; i<=endPos; i++) {
					Mutation mutation = new Mutation();
					String alt = alts.charAt(idx)+"";
					String ref = refs.charAt(idx)+"";
					mutation.altSeq = alt; mutation.refSeq = ref; mutation.pos = i; mutation.chr = chr;
					mutation.type = InputConvertorConstants.SNP;
					mutation.key = key;
					this.putMutation(chr, mutation);
					idx++;
				}
				totalOfMutations[InputConvertorConstants.SNP]++;
			} 
			// Deletion
			else if(refLen > altLen) {
				int idx = 0;
				for(int i=startPos+1; i<=endPos; i++) {
					Mutation mutation = new Mutation();
					String alt = InputConvertorConstants.DELETION_MARK;
					String ref = refs.charAt(idx)+"";
					mutation.altSeq = alt; mutation.refSeq = ref; mutation.pos = i; mutation.chr = chr;
					mutation.type = InputConvertorConstants.DEL;
					mutation.key = key;
					this.putMutation(chr, mutation);
					idx++;
				}
				totalOfMutations[InputConvertorConstants.DEL]++;
			} 
			// Insertion
			else if(refLen < altLen) {
				Mutation mutation = new Mutation();
				String alt = alts;
				String ref = InputConvertorConstants.DELETION_MARK;
				mutation.altSeq = alt; mutation.refSeq = ref; mutation.pos = startPos; mutation.chr = chr;
				mutation.type = InputConvertorConstants.DEL;
				mutation.key = key;
				this.putMutation(chr, mutation);
				totalOfMutations[InputConvertorConstants.INS]++;
			}
			
		}
		
		BR.close();
		
		System.out.println("A total of mutations: "+ (totalOfMutations[InputConvertorConstants.SNP] + totalOfMutations[InputConvertorConstants.INS] + totalOfMutations[InputConvertorConstants.DEL]));
		System.out.println("SNPs: "+totalOfMutations[InputConvertorConstants.SNP]);
		System.out.println("INSs: "+totalOfMutations[InputConvertorConstants.INS]);
		System.out.println("DELs: "+totalOfMutations[InputConvertorConstants.DEL]);
	}
	
	private void putMutation (String chr, Mutation mutation) {
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
	}
	
	
	private ArrayList<Mutation> getMutationByCondition (String chr, int start, int end, byte typeFlag) {
		ArrayList<Mutation> mutations = new ArrayList<Mutation>();
		if(positionalMap.get(chr) != null) {
			SortedMap<Integer, ArrayList<Mutation>> sortedMap = positionalMap.get(chr).subMap(start, end);
			if(sortedMap != null) {
				final byte flag = typeFlag;
				sortedMap.forEach((pos, m)->{
					if(flag == InputConvertorConstants.ALL_MUT) {
						mutations.addAll(m);
					} else {
						for(Mutation mutation :m) {
							if(mutation.type == flag) {
								mutations.add(mutation);
							}
						}
					}
					
				});
			}
		}
		
		
		return mutations;
	}
	
	
	/**
	 * zero-based <br>
	 * [start, end)
	 * @param start
	 * @param end
	 * @return
	 */
	public ArrayList<Mutation> getMutationByRange (String chr, int start, int end) {
		return getMutationByCondition(chr, start, end, InputConvertorConstants.ALL_MUT);
	}
	public ArrayList<Mutation> getSNPByRange (String chr, int start, int end) {
		return getMutationByCondition(chr, start, end, InputConvertorConstants.SNP);
	}
	public ArrayList<Mutation> getINSByRange (String chr, int start, int end) {
		return getMutationByCondition(chr, start, end, InputConvertorConstants.INS);
	}
	public ArrayList<Mutation> getDELByRange (String chr, int start, int end) {
		return getMutationByCondition(chr, start, end, InputConvertorConstants.DEL);
	}
}
