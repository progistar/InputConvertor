package zhanglab.inputconvertor.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

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
	
	public void setSequence(String chr, Exon[] exonGraph) {
		StringBuilder sequence = genomeMap.get(chr);
		ArrayList<Mutation> snps = new ArrayList<Mutation>();
		ArrayList<Mutation> inss = new ArrayList<Mutation>();
		ArrayList<Mutation> dels = new ArrayList<Mutation>();
		Exon startExon = exonGraph[0];
		Exon endExon = exonGraph[1];
		Exon nextExon = null;
		
		nextExon = startExon;
		
		while(nextExon != null) {
			// end of exon
			if(nextExon.start == Integer.MAX_VALUE) {
				nextExon = null;
			} 
			// exons in middle
			else if (nextExon.start != -1) {
				int start = nextExon.start - 1;
				int end = nextExon.end;
				String refSequence = sequence.subSequence(start,  end).toString();
				nextExon.refNucleotide = refSequence;
				if(vars != null) {
					snps.addAll(vars.getSNPByRange(chr, nextExon.start, nextExon.end+1));
					snps.addAll(vars.getMNPByRange(chr, nextExon.start, nextExon.end+1));
					
					dels.addAll(vars.getDELByRange(chr, nextExon.start, nextExon.end+1));
					inss.addAll(vars.getINSByRange(chr, nextExon.start, nextExon.end+1));
				}
				
				nextExon = nextExon.nextExons.get(0);
			}
			// start exon
			else {
				nextExon = nextExon.nextExons.get(0);
			}
		}
		
		// deletions
		for(Mutation del :dels) {
			nextExon = startExon;
			
			int delStartPos = del.pos;
			int delEndPos = del.pos + del.refSeq.length()-1;
			
			// connect delStartExon -> delExon -> delEndExon
			Exon delExon = new Exon(chr, delStartPos, delEndPos);
			delExon.type = del.type;
			delExon.refNucleotide = del.refSeq;
			delExon.altNucleotide = del.altSeq;
			
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
				
				nextExon = Exon.getWildTypeExon(nextExon.nextExons);
			}
			
			if(delEndExon == null) {
				delEndExon = endExon;
			}
			
			delStartExon.nextExons.add(delExon);
			delExon.prevExons.add(delStartExon);
			delEndExon.prevExons.add(delExon);
			delExon.nextExons.add(delEndExon);
		}
		
		// insertions
		for(Mutation ins :inss) {
			nextExon = startExon;
			int insStartPos = ins.pos;
			int insEndPos = ins.pos;
			
			// connect delStartExon -> delExon -> delEndExon
			Exon insExon = new Exon(chr, insStartPos, insEndPos);
			insExon.type = ins.type;
			insExon.refNucleotide = ins.refSeq;
			insExon.altNucleotide = ins.altSeq;
			
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
				
				nextExon = Exon.getWildTypeExon(nextExon.nextExons);
			}
			
			if(insEndExon == null) {
				insEndExon = endExon;
			}
			
			insStartExon.nextExons.add(insExon);
			insExon.prevExons.add(insStartExon);
			insEndExon.prevExons.add(insExon);
			insExon.nextExons.add(insEndExon);
		}
		
		// snps
		for(Mutation snp :snps) {
			nextExon = startExon;
			
			int snpStartPos = snp.pos;
			int snpEndPos = snp.pos + snp.refSeq.length()-1;
			
			Exon snpStartExon = null;
			Exon snpEndExon = null;
			
			while((nextExon.start != Integer.MAX_VALUE)) {
				if(nextExon.start < snpStartPos) {
					snpStartExon = nextExon;
				}
				if(nextExon.start > snpEndPos && snpEndExon == null) {
					snpEndExon = nextExon;
				}
				
				// if the exon has the position (mutation)
				if(nextExon.hasPosition(snpStartPos-1)) {
					nextExon = Exon.divdeExon(nextExon, snpStartPos-1);
					snpStartExon = nextExon;
				}
				if(nextExon.hasPosition(snpEndPos)) {
					nextExon = Exon.divdeExon(nextExon, snpEndPos);
				}
				
				nextExon = Exon.getWildTypeExon(nextExon.nextExons);
			}
			
			
			// coordinate snpExon by wildExon
			snpStartExon = Exon.getWildTypeExon(snpStartExon.nextExons);
			snpEndExon = Exon.getWildTypeExon(snpEndExon.prevExons);
			
			nextExon = snpStartExon;
			
			ArrayList<Exon> mnps = new ArrayList<Exon>();
			
			while(nextExon != null) {

				// connect delStartExon -> delExon -> delEndExon
				Exon snpExon = new Exon(chr, snpStartPos, snpEndPos);
				snpExon.type = snp.type;
				
				int shiftStartPos = nextExon.start - snpExon.start;
				int shiftEndPos = shiftStartPos + nextExon.refNucleotide.length();
				
				snpExon.start = nextExon.start; snpExon.end = nextExon.end;
				snpExon.refNucleotide = nextExon.refNucleotide;
				snpExon.altNucleotide = snp.altSeq.substring(shiftStartPos, shiftEndPos);
				mnps.add(snpExon);
				
				if(nextExon == snpEndExon) {
					nextExon = null;
				} else {
					nextExon = Exon.getWildTypeExon(nextExon.nextExons);
				}
			}
			
			ArrayList<Exon> mExons = Exon.mergeAdjacentExons(mnps);
			Exon firstExon = mExons.get(0);
			Exon lastExon = mExons.get(mExons.size()-1);
			
			for(Exon pExon : snpStartExon.prevExons) {
				pExon.nextExons.add(firstExon);
				firstExon.prevExons.add(pExon);
			}
			
			for(Exon nExon : snpEndExon.nextExons) {
				nExon.prevExons.add(lastExon);
				lastExon.nextExons.add(nExon);
			}
		}
		
		
	}
}
