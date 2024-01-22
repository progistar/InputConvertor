package zhanglab.inputconvertor;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.junit.jupiter.api.Test;

import junit.framework.Assert;
import zhanglab.inputconvertor.data.FastaEntry;
import zhanglab.inputconvertor.data.GTFLoader;
import zhanglab.inputconvertor.data.GenomeLoader;
import zhanglab.inputconvertor.input.CIRIquant;

public class CIRIquantTest {
	
	public File genomeFile = null;
	public File referenceFile = null;
	public File ciriquantFile = null;
	public ArrayList<String> expectedTranslations = null;
	//File vepFile = new File("/Users/seunghyukchoi/Documents/1_Projects/2023_Neoflow2/2_iRNAseq/vep/C3L-01632.txt");
	
	public CIRIquantTest () {
		this.loadTestData();
	}
	
	@Test
	public void translation () throws IOException {
		
		CIRIquantTest test = new CIRIquantTest();
		
		CIRIquant ciriquant = new CIRIquant(test.ciriquantFile);
		GTFLoader reference = new GTFLoader(test.referenceFile);
		GenomeLoader gLoader = new GenomeLoader(test.genomeFile);
		
		System.out.println(reference.geneToTranscripts.size());
		
		ciriquant.enrollGenomeSequence(gLoader);
		ciriquant.enrollReferenceGTF(reference);
		
		ArrayList<FastaEntry> entries = ciriquant.getFastaEntry();
		System.out.println(entries.size());
		
		BufferedWriter BW = new BufferedWriter(new FileWriter("test/test.out"));
		
        for(FastaEntry entry : entries) {
        	BW.append(">"+entry.toHeader());
        	BW.newLine();
        	BW.append(entry.sequence);
        	BW.newLine();
        }
        
        for(int i=0; i<24; i++) {
        	FastaEntry entry = entries.get(i);
        	Assert.assertTrue("Translated sequence was not equal: "+i + test.expectedTranslations.get(i) + "=>" + entry.sequence
        			,test.expectedTranslations.get(i).equalsIgnoreCase(entry.sequence));
        }
        
        BW.close();
	}

	public void loadTestData () {
		this.genomeFile = new File("test/test.genome.fa");
		this.referenceFile = new File("test/ref.gtf");
		this.ciriquantFile = new File("test/CIRIquant.gtf");
		this.expectedTranslations = new ArrayList<String>();
		
		// chr1:50|400_junction_exact
		this.expectedTranslations.add("PVCGDARHALPQHQVICCLLAQTSRVLST");
		this.expectedTranslations.add("QCVVMPGMPFPSIRLFAVSXPRLPVSFPP");
		this.expectedTranslations.add("SVWXCQACPSPASGYLLSLSPDFPCPFH");
		
		// chr1:200|400_junction_exact
		this.expectedTranslations.add("PVCGDARHALPQHQVSFPVERSHAXSGMG");
		this.expectedTranslations.add("QCVVMPGMPFPSIRXVFLWRGAMPRVGWA");
		this.expectedTranslations.add("SVWXCQACPSPASGEFSCGEEPCLEWDG");
		
		// chr1:200|250_junction_exact
		this.expectedTranslations.add("PVERSHAXSGMGHCSSFPVERSHAXSGMG");
		this.expectedTranslations.add("LWRGAMPRVGWAIVRVFLWRGAMPRVGWA");
		this.expectedTranslations.add("CGEEPCLEWDGPLFEFSCGEEPCLEWDG");
		
		// chr1:100|400_junction_shift
		this.expectedTranslations.add("PVCGDARHALPQHQVERSQGLDAVVFICR");
		this.expectedTranslations.add("QCVVMPGMPFPSIRLRGHRVLMLWSSSAG");
		this.expectedTranslations.add("SVWXCQACPSPASGXEVTGSXCCGLHLQ");
		
		// chr1:100|350_junction_shift
		this.expectedTranslations.add("KIGGKMSESINFSHNERSQGLDAVVFICR");
		this.expectedTranslations.add("RLEERXVRASTSLTMRGHRVLMLWSSSAG");
		this.expectedTranslations.add("DWRKDEXEHQLLSQXEVTGSXCCGLHLQ");
		
		// chr1:200|350_junction_shift
		this.expectedTranslations.add("KIGGKMSESINFSHKSFPVERSHAXSGMG");
		this.expectedTranslations.add("RLEERXVRASTSLTRVFLWRGAMPRVGWA");
		this.expectedTranslations.add("DWRKDEXEHQLLSQEFSCGEEPCLEWDG");

		// chr1:280|400_intron
		this.expectedTranslations.add("PVCGDARHALPQHQVLIPQPGIGERLEER");
		this.expectedTranslations.add("QCVVMPGMPFPSIRSXYHNQAXGKDWRKD");
		this.expectedTranslations.add("SVWXCQACPSPASGLNTTTRHRGKIGGK");
		
		// chr1:500|600_junction_exact
		this.expectedTranslations.add("GLDLSPGGGQSHLWFSIVLLDQXYTRHPV");
		this.expectedTranslations.add("AWIXALVEVKATFGSAXCSWTSDTPGTLS");
		this.expectedTranslations.add("PGSEPWWRSKPPLVQHSAPGPVIHPAPC");
		
		// chr1:550|580_junction_shift
		this.expectedTranslations.add("EPWWTLLAWIXALVDAVGLDLSPGGRCWP");
		this.expectedTranslations.add("SPGGRCWPGSEPWWTLLAWIXALVDAVGL");
		this.expectedTranslations.add("ALVDAVGLDLSPGGRCWPGSEPWWTLLA");
		
		
		
	}
	
}
