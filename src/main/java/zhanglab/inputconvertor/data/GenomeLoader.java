package zhanglab.inputconvertor.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Hashtable;
import java.util.LinkedList;

import zhanglab.inputconvertor.env.InputConvertorConstants;

public class GenomeLoader {

	private Hashtable<String, StringBuilder> genomeMap = new Hashtable<String, StringBuilder>();
	private VEPLoader vars = null;
	
	public void clear() {
		genomeMap.clear();
	}
	
	public GenomeLoader (File file) {
		// load genome
		try {
			
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String line = null;
			
			StringBuilder sequence = null;
			while((line = BR.readLine()) != null) {
				if(line.startsWith(">")) {
					String chr = line.split("\\s")[0].substring(1);
					System.out.print("Read "+chr+"... ");
					sequence = new StringBuilder();
					genomeMap.put(chr, sequence);
				} else {
					sequence.append(line.toUpperCase());
				}
			}
			
			BR.close();
			System.out.println();
		}catch(IOException ioe) {
			
		}
	}
	
	public void enrollVEPLaoder (VEPLoader vep) {
		System.out.println("## GenomeLoader ##");
		System.out.println("Enroll Somatic VEP");
		this.vars = vep;
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
	
	public void setSequence(String chr, ArrayList<Exon> exons) {
		StringBuilder sequence = genomeMap.get(chr);
		ArrayList<Mutation> snps = new ArrayList<Mutation>();
		ArrayList<Mutation> inss = new ArrayList<Mutation>();
		ArrayList<Mutation> dels = new ArrayList<Mutation>();
		for(Exon exon : exons) {
			// start and end exons
			if(exon.start == -1 || exon.start == Integer.MAX_VALUE) continue;
			int start = exon.start - 1;
			int end = exon.end;
			String refSequence = sequence.subSequence(start,  end).toString();
			exon.nucleotide = refSequence;
			if(vars != null) {
				snps.addAll(vars.getSNPByRange(chr, exon.start, exon.end+1));
				dels.addAll(vars.getDELByRange(chr, exon.start, exon.end+1));
				inss.addAll(vars.getINSByRange(chr, exon.start, exon.end+1));
			}
		}
		
		// deletions
		for(Mutation del :dels) {
			Exon nextExon = exons.get(0);
			
			int delStartPos = del.pos;
			int delEndPos = del.pos + del.refSeq.length()-1;
			
			// connect delStartExon -> delExon -> delEndExon
			Exon delExon = new Exon(chr, delStartPos, delEndPos);
			delExon.mutation = del;
			delExon.type = InputConvertorConstants.DEL;
			
			Exon delStartExon = null;
			Exon delEndExon = null;
			while((nextExon.start != Integer.MAX_VALUE)) {
				if(nextExon.start < delStartPos) {
					delStartExon = nextExon;
				}
				if(nextExon.start > delEndPos && delEndExon == null) {
					delEndExon = nextExon;
				}
				
				// if the exon has the position (mutation)
				if(nextExon.hasPosition(delStartPos-1)) {
					nextExon = Exon.divdeExon(nextExon, delStartPos-1);
					delStartExon = nextExon;
				}
				if(nextExon.hasPosition(delEndPos)) {
					nextExon = Exon.divdeExon(nextExon, delEndPos);
				}
				
				Exon refExon = null;
				
				for(Exon nExon : nextExon.nextExons) {
					if(nExon.mutation == null) {
						refExon = nExon;
						break;
					}
				}
				nextExon = refExon;
			}
			
			if(delEndExon == null) {
				delEndExon = exons.get(exons.size()-1);
			}
			
			delStartExon.nextExons.add(delExon);
			delExon.prevExons.add(delStartExon);
			delEndExon.prevExons.add(delExon);
			delExon.nextExons.add(delEndExon);
		}
		
		// insertions
		for(Mutation ins :inss) {
			Exon nextExon = exons.get(0);
			int insStartPos = ins.pos;
			int insEndPos = ins.pos;
			
			// connect delStartExon -> delExon -> delEndExon
			Exon insExon = new Exon(chr, insStartPos, insEndPos);
			insExon.mutation = ins;
			insExon.type = InputConvertorConstants.INS;
			
			Exon insStartExon = nextExon;
			Exon insEndExon = null;
			while((nextExon.start != Integer.MAX_VALUE)) {
				if(nextExon.start > insEndPos && insEndExon == null) {
					insEndExon = nextExon;
				}
				// if the exon has the position (mutation)
				if(nextExon.hasPosition(insStartPos)) {
					nextExon = Exon.divdeExon(nextExon, insStartPos);
					insStartExon = nextExon;
				}
				
				Exon refExon = null;
				
				for(Exon nExon : nextExon.nextExons) {
					if(nExon.type == InputConvertorConstants.WILD) {
						refExon = nExon;
						break;
					}
				}
				nextExon = refExon;
			}
			
			if(insEndExon == null) {
				insEndExon = exons.get(exons.size()-1);
			}
			
			insStartExon.nextExons.add(insExon);
			insExon.prevExons.add(insStartExon);
			insEndExon.prevExons.add(insExon);
			insExon.nextExons.add(insEndExon);
		}
		
		// snps
		for(Mutation snp :snps) {
			Exon nextExon = exons.get(0);
			
			int snpStartPos = snp.pos;
			int snpEndPos = snp.pos;
			
			// connect delStartExon -> delExon -> delEndExon
			Exon snpExon = new Exon(chr, snpStartPos, snpEndPos);
			snpExon.mutation = snp;
			snpExon.type = InputConvertorConstants.SNP;
			
			while((nextExon.start != Integer.MAX_VALUE)) {
				// if the exon has the position (mutation)
				if(nextExon.hasPosition(snpStartPos)) {
					nextExon = Exon.divdeExon(nextExon, snpStartPos);
					
					if(nextExon.hasPosition(snpStartPos-1)) {
						nextExon = Exon.divdeExon(nextExon, snpStartPos-1);
						nextExon = nextExon.nextExons.get(0);
					}
					
					for(Exon pExon : nextExon.prevExons) {
						pExon.nextExons.add(snpExon);
						snpExon.prevExons.add(pExon);
					}
					
					for(Exon nExon : nextExon.nextExons) {
						nExon.prevExons.add(snpExon);
						snpExon.nextExons.add(nExon);
					}
					
					// found
					break;
				}
				
				for(Exon nExon : nextExon.nextExons) {
					if(nExon.type == InputConvertorConstants.WILD) {
						nextExon = nExon;
						break;
					}
				}
			}
		}
		
		
	}
}
