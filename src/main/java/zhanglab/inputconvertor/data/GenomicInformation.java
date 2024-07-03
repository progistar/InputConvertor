package zhanglab.inputconvertor.data;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;

import zhanglab.inputconvertor.env.InputConvertorConstants;

public class GenomicInformation {
	
	public String fastaId;
	public String mark;
	public String ensgId;
	public String enstId;
	public String geneName;
	public String strand;
	public String location;
	public String ilPeptide; // distinguish IL
	public String from;
	public int frame;
	public int penalty;
	public int mutCnt;
	public boolean isContam;;
	
	public GenomicInformation (String header, Emit emit, String fullSequence, boolean isContam) {
		this.isContam = isContam;
		
		if(isContam) {
			this.fastaId = header.substring(1);
			this.mark = InputConvertorConstants.CONT_HEADER_ID;
			this.ensgId = "Contaminant";
			this.enstId = "Contaminant";
			this.geneName = "Contaminant";
			this.frame = 0;
			this.strand = "+";
			this.ilPeptide = fullSequence.substring(emit.getStart(), emit.getEnd()+1);
			this.mutCnt = 0;
			this.location = "-";
			this.penalty = 0;
			this.from = "Contaminant";
			return;
		}
		
		String[] fields = header.split("\\|");
		this.fastaId = fields[InputConvertorConstants.FASTA_ID_IDX].substring(1);
		this.mark = this.fastaId.split("\\_")[0].replaceAll("[0-9]", "");
		this.ensgId = fields[InputConvertorConstants.ENSG_ID_IDX];
		this.enstId = fields[InputConvertorConstants.ENST_ID_IDX];
		this.geneName = fields[InputConvertorConstants.GENE_NAME_IDX];
		this.strand = fields[InputConvertorConstants.STRAND_IDX];
		this.ilPeptide = fullSequence.substring(emit.getStart(), emit.getEnd()+1);
		this.mutCnt = 0;
		
		
		
		// calculate genomic Loci
		// if the sequence comes from fusion genes or circRNAs then it uses the reported genomic location
		if(this.mark.startsWith(InputConvertorConstants.ARRIBA_HEADER_ID) ||
			this.mark.startsWith(InputConvertorConstants.CIRIQUANT_HEADER_ID)) {
			this.location = fields[InputConvertorConstants.LOCATION_IDX];
		} else {
			this.frame = Integer.parseInt(fields[InputConvertorConstants.FRAME_IDX]);
			String[] genomicLoci = fields[InputConvertorConstants.LOCATION_IDX].split("@");
			int rnaLen = 0;
			for(int i=0; i<genomicLoci.length; i++) {
				String[] varInfo = genomicLoci[i].split("\\[");
				String[] gInfo = varInfo[0].split("\\:");
				// [start, end] one-based
				int start = Integer.parseInt(gInfo[1].split("\\-")[0]);
				int end = Integer.parseInt(gInfo[1].split("\\-")[1]);
				
				rnaLen += (end-start+1);
			}
			
			
			// [ startSize, endSize] zero-based
			int startSize = 0;
			int endSize = 0;
			
			startSize = emit.getStart();
			endSize = (emit.getEnd()+1);
			
			startSize = startSize*3;
			endSize = endSize*3-1;
			
			startSize += this.frame;
			endSize += this.frame;
			
			// forward strand
			if(strand.equalsIgnoreCase("+")) {
				
			} 
			// reverse strand
			else {
				startSize = rnaLen - startSize - 1;
				endSize = rnaLen - endSize - 1;
				int tmp = startSize;
				startSize = endSize;
				endSize = tmp;
			}
			
			

			ArrayList<Integer> starts = new ArrayList<Integer>();
			ArrayList<Integer> ends = new ArrayList<Integer>();
			ArrayList<String> vars = new ArrayList<String>();
			
			String chr = null;
			for(int i=0; i<genomicLoci.length; i++) {
				String[] varInfo = genomicLoci[i].split("\\[");
				String var = "";
				if(varInfo.length > 1) {
					var = "["+varInfo[1];
				}
				
				String[] gInfo = varInfo[0].split("\\:");
				chr = gInfo[0];
				// [start, end] one-based
				int start = Integer.parseInt(gInfo[1].split("\\-")[0]);
				int end = Integer.parseInt(gInfo[1].split("\\-")[1]);
				
				if(start+startSize > end) {
					startSize -= (end-start+1);
					endSize -= (end-start+1);
				} else {
					vars.add(var);
					if(startSize == -1) {
						starts.add(start);
					} else {
						starts.add(start+startSize);
						startSize = -1;
					}
					
					if(start+endSize > end) {
						ends.add(end);
						endSize -= (end-start+1);
					} else {
						ends.add(start + endSize);
						endSize = -1;
						break;
					}
				}
			}
			StringBuilder str = new StringBuilder();
			for(int i=0; i<starts.size(); i++) {
				if(i!=0) {
					str.append("|");
				}
				str.append(chr).append(":").append(starts.get(i)).append("-").append(ends.get(i)).append(vars.get(i));
			}
			this.location = str.toString();
		}
		
		// scoring

		// mutation starts with '[' 
		for(int i=0; i<this.location.length(); i++) {
			if(this.location.charAt(i) == '[') {
				this.mutCnt ++;
			}
		}
		
		// penalty
		if(this.mark.startsWith(InputConvertorConstants.ARRIBA_HEADER_ID) ||
				this.mark.startsWith(InputConvertorConstants.CIRIQUANT_HEADER_ID)) {
			this.penalty = 100000;
			
			if(this.mark.startsWith(InputConvertorConstants.ARRIBA_HEADER_ID)) {
				this.from = "Fusion gene";
			}
			if(this.mark.startsWith(InputConvertorConstants.CIRIQUANT_HEADER_ID)) {
				this.from = "CircRNA";
			}
			
		} else if(this.mark.startsWith(InputConvertorConstants.IRFINDER_HEADER_ID)) {
			this.penalty = 10000;
			this.from = "IR";
		} else if(this.mark.startsWith(InputConvertorConstants.STRINGTIE_HEADER_ID)) {
			this.penalty = 1000;
			this.from = "StringTie";
		} else if(this.mark.startsWith(InputConvertorConstants.EXON_TRANSLATION_HEADER_ID)) {
			this.penalty = 100;
			this.from = "Exon";
		} else if(this.mark.startsWith(InputConvertorConstants.REF_HEADER_ID)) {
			this.penalty = 10;
			this.from = "Reference";
		}
		this.penalty += this.mutCnt;
		
	}
	
	public String getKey () {
		return this.location;
	}
	
	
	
	public static void main(String[] args) throws IOException {
		String header = ">@RF26660_DB_C3L-00973.T|ENSG00000157856.12|ENST00000288710.7|DRC1|0|+|chr2:26401990-26402144@chr2:26414344-26414431@chr2:26421288-26421400@chr2:26424271-26424454@chr2:26429628-26429765@chr2:26430786-26430872@chr2:26431884-26432006@chr2:26440378-26440517@chr2:26444222-26444356@chr2:26444716-26444948@chr2:26448691-26448803@chr2:26449996-26450085@chr2:26450592-26450681@chr2:26453320-26453549@chr2:26454647-26454790@chr2:26455131-26455233@chr2:26456461-26456514";
		String sequence = "MNPPGSLEALDPNVDEHLSTQILAPSVHSDNSQERIQARRLRIAARLEARRREALGEYLDGKKESEEDQSKSYKQKEESRLKLAKLLLCGTELVTNIQVAIDIREIHRRVEEEEIKRQRIEKLENEVKTSQDKFDEITSKWEEGKQKRIPQELWEMLNTQQLHCAGLLEDKNKLISELQQELKTKDDQYVKDLKKQSDDICLLLERMEEQVKNVMKTFREELYNIEKAFEVERQELLASNKKKWEQALQAHNAKELEYLNNRMKKVEDYEKQLNRQRIWDCEEYNMIKIKLEQDVQILEQQLQQRKAIYQLNQEKLEYNLQVLKKRDEESTVIKSQQKRKINRLHDILNNLRSKYAKQIKQFQEENQSLTSDYKRLVMQFKELQKAMRHFALIDDEKFWEIWLMNEEEAKDLIARAFDVDRIIHTHHLGLPWAAPDFWFLNNVGPISQQPQKSATQIVEEMLMRSEEEEAEEAAAEPESYLDLPKQISEKTTKRILMLLCDESGFLIESKLLSLLLPLEQNECYLLRLDAIFSALGIESEDDLYKLVNFFLKYRAHRLSSSLQIKPCSQASMEKASMEETSTRSELELAEQTEMEGEKEESLVEGEKEEEEETPPSPWVIHPNDVLKILEAFVMGLKKPRDSRAPLRVQKNVRDNSKDSEYWQALTTVIPSSKQNLWDALYTALEKYHLVLTQRAKLLLENSSLEQQNTELQALLQQYLNSKINSELQVPPTQVLRVPTK";
		
		String headerRV = ">@RF41428_DB_C3L-00973.T|ENSG00000117139.18|ENST00000648473.1|KDM5B|0|-|chr1:202729457-202729474@chr1:202729707-202730027@chr1:202730909-202731063@chr1:202731828-202731939@chr1:202733401-202733886@chr1:202735429-202735587@chr1:202736213-202736392@chr1:202740674-202740812@chr1:202741367-202741722@chr1:202742391-202742505@chr1:202742655-202742805@chr1:202745858-202745982@chr1:202746142-202746323@chr1:202748945-202749139@chr1:202750659-202750778@chr1:202752905-202753067@chr1:202755271-202755452@chr1:202756358-202756516@chr1:202758391-202758510@chr1:202760415-202760573@chr1:202762699-202762808@chr1:202764049-202764145@chr1:202766926-202767060@chr1:202773118-202773288@chr1:202774613-202774735@chr1:202777017-202777094@chr1:202808102-202808305";
		String sequenceRV = "MEAATTLHPGPRPALPLGGPGPLGEFLPPPECPVFEPSWEEFADPFAFIHKIRPIAEQTGICKVRPPPDWQPPFACDVDKLHFTPRIQRLNELEAQTRVKLNFLDQIAKYWELQGSTLKIPHVERKILDLFQLNKLVAEEGGFAVVCKDRKWTKIATKMGFAPGKAVGSHIRGHYERILNPYNLFLSGDSLRCLQKPNLTTDTKDKEYKPHDIPQRQSVQPSETCPPARRAKRMRAEAMNIKIEPEETTEARTHNLRRRMGCPTPKCENEKEMKSSIKQEPIERKDYIVENEKEKPKSRSKKATNAVDLYVCLLCGSGNDEDRLLLCDGCDDSYHTFCLIPPLHDVPKGDWRCPKCLAQECSKPQEAFGFEQAARDYTLRTFGEMADAFKSDYFNMPVHMVPTELVEKEFWRLVSTIEEDVTVEYGADIASKEFGSGFPVRDGKIKLSPEEEEYLDSGWNLNNMPVMEQSVLAHITADICGMKLPWLYVGMCFSSFCWHIEDHWSYSINYLHWGEPKTWYGVPGYAAEQLENVMKKLAPELFVSQPDLLHQLVTIMNPNTLMTHEVPVYRTNQCAGEFVITFPRAYHSGFNQGFNFAEAVNFCTVDWLPLGRQCVEHYRLLHRYCVFSHDEMICKMASKADVLDVVVASTVQKDMAIMIEDEKALRETVRKLGVIDSERMDFELLPDDERQCVKCKTTCFMSAISCSCKPGLLVCLHHVKELCSCPPYKYKLRYRYTLDDLYPMMNALKLRAESYNEWALNVNEALEAKINKKKSLVSFKALIEESEMKKFPDNDLLRHLRLVTQDAEKCASVAQQLLNGKRQTRYRSGGGKSQNQLTVNELRQFVTQLYALPCVLSQTPLLKDLLNRVEDFQQHSQKLLSEETPSAAELQDLLDVSFEFDVELPQLAEMRIRLEQARWLEEVQQACLDPSSLTLDDMRRLIDLGVGLAPYSAVEKAMARLQELLTVSEHWDDKAKSLLKARPRHSLNSLATAVKEIEEIPAYLPNGAALKDSVQRARDWLQDVEGLQAGGRVPVLDTLIELVTRGRSIPVHLNSLPRLETLVAEVQAWKECAVNTFLTENSPYSLLEVLCPRCDIGLLGLKRKQRKLKEPLPNGKKKSTKLESLSDLERALTESKETASAMATLGEARLREMEALQSLRLANEGKLLSPLQDVDIKICLCQKAPAAPMIQCELCRDAFHTSCVAVPSISQGLRIWLCPHCRRSEKPPLEKILPLLASLQRIRVRLPEGDALRYMIERTVNWQHRAQQLLSSGNLKFVQDRVGSGLLYSRWQASAGQVSDTNKVSQPPGTTSFSLPDDWDNRTSYLHSPFSTGRSCIPLHGVSPEVNELLMEAQLLQVSLPEIQELYQTLLAKPSPAQQTDRSSPVRPSSEKNDCCRGKRDGINSLERKLKRRLEREGLSSERWERVKKMRTPKKKKIKLSHPKDMNNFKLERERSYELVRSAETHSLPSDTSYSEQEDSEDEDAICPAVSCLQPEGDELSPSRR";
		
		String peptide1 = "MNPPGSLEAL";
		String peptide2 = "RLEARRREALGEYLDGKKESEE";
		String peptide3 = "TNIQVAIDIREIHRRVEEEEIKRQRIEKLENEVKTSQDKFDEITSKWEEGKQKRIPQELWEMLNTQQLHCAGLLEDKNKLISELQQELKTKDDQ";
		
		String peptide4 = "MEAATTLHPGPRPALP";
		String peptide5 = "LPPPECPVFEPSWEEFADPFAFIHKIRPIAEQTGICKVRPPPDWQPPFACDVDKL";
		String peptide6 = "CKDRKWTKIATKMGFAPGKAVGSHIRGHYERILNPYNLFLSGDSLRCLQKPNLTTDTKDKEYKPH";
		
		String headerVar = ">@EX@VR2710_DB_C3L-00973.T|ENSG00000040731.10|ENST00000264463.8|CDH10|0|-|chr5:24492817-24492881@chr5:24492882-24492882[SNP]C>A@chr5:24492883-24492925";
		String sequenceVar = "LIQTISAVDKDDPLVGQKFFFSLAAVNPNFTVQDNE";
		
		String peptide7 = "LIQTISAVDKDDPLVGQ";
		String peptide8 = "FFSLAAVNPNFTVQDNE";
		
		
		String headerVar2 = ">@EX@VR2711_DB_C3L-00973.T|ENSG00000040731.10|ENST00000264463.8|CDH10|1|-|chr5:24492817-24492881@chr5:24492882-24492882[SNP]C>A@chr5:24492883-24492925";
		String sequenceVar2 = "XYRLXVQXTKMTLXLDRNFFSVXLLSIQTSQYRIMK";
		
		String peptide9 = "LDRNFFSV";
		
		Trie trie = Trie.builder().addKeyword(peptide1).addKeyword(peptide2).addKeyword(peptide3)
				.addKeyword(peptide4).addKeyword(peptide5).addKeyword(peptide6)
				.addKeyword(peptide7).addKeyword(peptide8)
				.addKeyword(peptide9).build();
		
		Collection<Emit> emits = trie.parseText(sequence);
		
		System.out.println(header);
		for(Emit emit : emits) {
			GenomicInformation a = new GenomicInformation(header, emit, sequence, false);
		}
		
		emits = trie.parseText(sequenceRV);
		System.out.println(headerRV);
		for(Emit emit : emits) {
			GenomicInformation a = new GenomicInformation(headerRV, emit, sequenceRV, false);
		}
		
		emits = trie.parseText(sequenceVar);
		System.out.println(headerVar);
		for(Emit emit : emits) {
			GenomicInformation a = new GenomicInformation(headerVar, emit, sequenceVar, false);
		}
		
		emits = trie.parseText(sequenceVar2);
		System.out.println(headerVar2);
		for(Emit emit : emits) {
			GenomicInformation a = new GenomicInformation(headerVar2, emit, sequenceVar2, false);
		}
	}
}
