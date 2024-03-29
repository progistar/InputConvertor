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
	
	public VEPLoader (File file) throws IOException {
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		String[] header = BR.readLine().split("\t");
		
		int locIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.VEP_LOC_FIELD_NAME);
		int refIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.VEP_REF_FIELD_NAME);
		int altIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.VEP_ALT_FIELD_NAME);
		
		int[] totalOfMutations = new int[10];
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
			
			String mark = ">";
			if(isSomatic) {
				mark = ">>";
			}
			byte type = InputConvertorConstants.WILD;
			String typeStr = null;
			if(refLen == altLen) {
				if(refLen == 1) {
					type = InputConvertorConstants.SNP;
					typeStr = "SNP";
				}
				// MNP (DNP , TNP ...)
				else {
					type = InputConvertorConstants.MNP;
					typeStr = "MNP";
				}
			} else if(refLen > altLen) {
				startPos++;
				type = InputConvertorConstants.DEL;
				typeStr = "DEL";
			} else if(refLen < altLen) {
				type = InputConvertorConstants.INS;
				typeStr = "INS";
			}
			
			String key = chr+":"+startPos+"-"+endPos+"["+typeStr+"]"+refs+mark+alts;
			if(removeDuplication.get(key) != null) {
				continue;
			}
			removeDuplication.put(key, "");
			
			//**
			// VEP can contain di/tri something like that
			
			// SNP / DNP / TNP ...
			if(type == InputConvertorConstants.SNP || type == InputConvertorConstants.MNP) {
				Mutation mutation = new Mutation();
				String alt = alts;
				String ref = refs;
				mutation.altSeq = alt; mutation.refSeq = ref; mutation.pos = startPos; mutation.chr = chr;
				mutation.type = type;
				mutation.key = key;
				this.putMutation(chr, mutation);
				totalOfMutations[mutation.type]++;
			} 
			// Deletion
			else if(type == InputConvertorConstants.DEL) {
				Mutation mutation = new Mutation();
				String alt = "";
				String ref = refs;
				mutation.altSeq = alt; mutation.refSeq = ref; mutation.pos = startPos; mutation.chr = chr;
				mutation.type = type;
				mutation.key = key;
				this.putMutation(chr, mutation);
				totalOfMutations[mutation.type]++;
			} 
			// Insertion
			else if(type == InputConvertorConstants.INS) {
				Mutation mutation = new Mutation();
				String alt = alts;
				String ref = "";
				mutation.altSeq = alt; mutation.refSeq = ref; mutation.pos = startPos; mutation.chr = chr;
				mutation.type = type;
				mutation.key = key;
				this.putMutation(chr, mutation);
				totalOfMutations[mutation.type]++;
			}
			
		}
		
		BR.close();
		
		System.out.println("A total of mutations: "+ 
				(totalOfMutations[InputConvertorConstants.SNP] +
				totalOfMutations[InputConvertorConstants.MNP] +
				totalOfMutations[InputConvertorConstants.INS] + 
				totalOfMutations[InputConvertorConstants.DEL]));
		System.out.println("SNPs: "+totalOfMutations[InputConvertorConstants.SNP]);
		System.out.println("MNPs: "+totalOfMutations[InputConvertorConstants.MNP]);
		System.out.println("INSs: "+totalOfMutations[InputConvertorConstants.INS]);
		System.out.println("DELs: "+totalOfMutations[InputConvertorConstants.DEL]);
	}
	
	private void putMutation (String chr, Mutation mutation) {
		TreeMap<Integer, ArrayList<Mutation>> map = positionalMap.get(chr);
		if(map == null) {
			map = new TreeMap<Integer, ArrayList<Mutation>>();
			positionalMap.put(chr, map);
		}
		
		int startPos = mutation.pos;
		int endPos = mutation.pos + mutation.refSeq.length() - 1;
		
		if(mutation.type == InputConvertorConstants.INS) {
			endPos = startPos;
		}
		
		for(int pos = startPos; pos <= endPos; pos++) {
			ArrayList<Mutation> mutations = map.get(pos);
			if(mutations == null) {
				mutations = new ArrayList<Mutation>();
				map.put(pos, mutations);
			}
			mutations.add(mutation);
		}
	}
	
	
	private ArrayList<Mutation> getMutationByCondition (String chr, int start, int end, byte typeFlag) {
		ArrayList<Mutation> mutations = new ArrayList<Mutation>();
		Hashtable<String, Boolean> isDup = new Hashtable<String, Boolean>();
		
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
								if(isDup.get(mutation.key) == null) {
									mutations.add(mutation);
									isDup.put(mutation.key, true);
								}
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
	public ArrayList<Mutation> getMNPByRange (String chr, int start, int end) {
		return getMutationByCondition(chr, start, end, InputConvertorConstants.MNP);
	}
	public ArrayList<Mutation> getINSByRange (String chr, int start, int end) {
		return getMutationByCondition(chr, start, end, InputConvertorConstants.INS);
	}
	public ArrayList<Mutation> getDELByRange (String chr, int start, int end) {
		return getMutationByCondition(chr, start, end, InputConvertorConstants.DEL);
	}
}
