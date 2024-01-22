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
import zhanglab.inputconvertor.input.IRFinder;

public class IRFinderTest {

	public File genomeFile = null;
	public File referenceFile = null;
	public File irfinderFile = null;
	public ArrayList<String> expectedTranslations = null;
	
	public IRFinderTest() {
		this.loadTestData();
	}
	
	@Test
	public void translation () throws IOException {
		
		IRFinderTest test = new IRFinderTest();
		
		IRFinder irfinder = new IRFinder(test.irfinderFile);
		GTFLoader reference = new GTFLoader(test.referenceFile);
		GenomeLoader gLoader = new GenomeLoader(test.genomeFile);
		
		irfinder.enrollGenomeSequence(gLoader);
		irfinder.enrollReferenceGTF(reference);
		
		ArrayList<FastaEntry> entries = irfinder.getFastaEntry();
		System.out.println(entries.size());
		
		BufferedWriter BW = new BufferedWriter(new FileWriter("test/test.out"));
		
        for(FastaEntry entry : entries) {
        	BW.append(">"+entry.toHeader());
        	BW.newLine();
        	BW.append(entry.sequence);
        	BW.newLine();
        }
        
        BW.close();
        
        
        for(int i=0; i<6; i++) {
        	FastaEntry entry = entries.get(i);
        	Assert.assertTrue("Translated sequence was not equal: "+i + test.expectedTranslations.get(i) + "=>" + entry.sequence
        			,test.expectedTranslations.get(i).equalsIgnoreCase(entry.sequence));
        }
        
        
	}

	public void loadTestData () {
		this.genomeFile = new File("test/test.genome.fa");
		this.referenceFile = new File("test/ref.gtf");
		this.irfinderFile = new File("test/IRFinder.tsv");
		this.expectedTranslations = new ArrayList<String>();
		
		// chr1:50|400_junction_exact
		this.expectedTranslations.add("PVERSHAXSGMGHCSSSGPCCLHVTXYHNQAXGKDWRKDEXEHQL");
		this.expectedTranslations.add("LWRGAMPRVGWAIVHLLAPVVCMXLNTTTRHRGKIGGKMSESINF");
		this.expectedTranslations.add("CGEEPCLEWDGPLFIFWPLLSACNLIPQPGIGERLEERXVRASTS");
		
		// chr1:200|400_junction_exact
		this.expectedTranslations.add("SSGQSKGRCSCLSGRGLLLVKLGRQKAVRNVISGW");
		this.expectedTranslations.add("PAGRAKEGAAACQEEAYFWXNWADKRQXEMXSRGG");
		this.expectedTranslations.add("QRAEQRKVQLPVRKRPTSGETGQTKGSEKCDLGVV");
	}
}
