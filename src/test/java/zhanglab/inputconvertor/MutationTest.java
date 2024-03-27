package zhanglab.inputconvertor;

import java.io.IOException;
import java.util.ArrayList;

import org.junit.jupiter.api.Test;

import zhanglab.inputconvertor.data.Exon;
import zhanglab.inputconvertor.data.Mutation;
import zhanglab.inputconvertor.data.Transcript;

public class MutationTest {

	public ArrayList<Exon> makeTestData() {
		Transcript transcript = new Transcript();
		ArrayList<Exon> exons = new ArrayList<Exon>();
		
		Exon e1 = new Exon("", 1, 10);
		e1.nucleotide = "AAAAAAAAAA";
		Exon e2 = new Exon("", 21, 30);
		e2.nucleotide = "AAAAAAAAAA";
		Exon e3 = new Exon("", 32, 41);
		e3.nucleotide = "AAAAAAAAAA";
		Exon e4 = new Exon("", 51, 60);
		e4.nucleotide = "AAAAAAAAAA";
		
		exons.add(e1); exons.add(e2); exons.add(e3); exons.add(e4);
		transcript.exons = exons;
		transcript.setExons();
		
		return Exon.addStartEndEonxs(transcript.exons);
		
	}
	
	@Test
	public void runTest () throws IOException {
		ArrayList<Mutation> dels = new ArrayList<Mutation>();
		Mutation del1 = new Mutation();
		del1.altSeq = "-";
		del1.refSeq = "AAT";
		del1.pos = 5;
		dels.add(del1);
		
		ArrayList<Mutation> inss = new ArrayList<Mutation>();
		Mutation ins1 = new Mutation();
		ins1.altSeq = "AAT";
		ins1.refSeq = "-";
		ins1.pos = 5;
		inss.add(ins1);
		
		ArrayList<Exon> originExons = makeTestData();

		// deletions
		for(Mutation del :dels) {
			Exon nextExon = originExons.get(0);
			
			int delStartPos = del.pos;
			int delEndPos = del.pos + del.refSeq.length()-1;
			
			// connect delStartExon -> delExon -> delEndExon
			Exon delExon = new Exon("", delStartPos, delEndPos);
			delExon.mutation = del;
			Exon delStartExon = null;
			Exon delEndExon = null;
			while((nextExon.start != Integer.MAX_VALUE)) {
				
				System.out.println(nextExon.start+"-"+nextExon.end);
				
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
				System.out.println(nextExon.start+"-"+nextExon.end);
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
				delEndExon = originExons.get(originExons.size()-1);
			}
			
			delStartExon.nextExons.add(delExon);
			delExon.prevExons.add(delStartExon);
			delEndExon.prevExons.add(delExon);
			delExon.nextExons.add(delEndExon);
		}
		
		

		// insertions
		for(Mutation ins :inss) {
			Exon nextExon = originExons.get(0);
			
			int insStartPos = ins.pos;
			int insEndPos = ins.pos;
			
			// connect delStartExon -> delExon -> delEndExon
			Exon insExon = new Exon("", insStartPos, insEndPos);
			insExon.mutation = ins;
			Exon insStartExon = null;
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
					if(nExon.mutation == null) {
						refExon = nExon;
						break;
					}
				}
				nextExon = refExon;
			}
			
			if(insEndExon == null) {
				insEndExon = originExons.get(originExons.size()-1);
			}
			
			insStartExon.nextExons.add(insExon);
			insExon.prevExons.add(insStartExon);
			insEndExon.prevExons.add(insExon);
			insExon.nextExons.add(insEndExon);
		}
		
		System.out.println("Done");
		
	}
}
