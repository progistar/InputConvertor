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
import zhanglab.inputconvertor.data.VEPLoader;
import zhanglab.inputconvertor.input.Reference;

public class MutationTest {

	public File genomeFile = null;
	public File referenceFile = null;
	public File vepFile = null;
	public ArrayList<String> expectedTranslations = null;
	
	@Test
	public void runTest () throws IOException {

		MutationTest test = new MutationTest();
		test.loadTestData();
		
		GTFLoader refGTF = new GTFLoader(test.referenceFile);
		GenomeLoader gLoader = new GenomeLoader(test.genomeFile);
		VEPLoader vepLoader = null;
		if(test.vepFile != null) {
			vepLoader = new VEPLoader(test.vepFile);
			gLoader.enrollVEPLaoder(vepLoader);
		}
		
		
		Reference reference = new Reference(refGTF);
        reference.enrollGenomeSequence(gLoader);
        ArrayList<FastaEntry> entries = reference.getFastaEntry(false);
		
		System.out.println(entries.size());
		
		BufferedWriter BW = new BufferedWriter(new FileWriter("test/test.out"));
		
        for(FastaEntry entry : entries) {
        	System.out.println(entry.description+"   " +entry.frame);
        	System.out.println(entry.sequence);
        	BW.append(">"+entry.toHeader("UNIQUE"));
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
		this.referenceFile = new File("test/comb_mut_ref.gtf");
		this.vepFile = new File("test/comb_mut_vep.tsv");
		this.expectedTranslations = new ArrayList<String>();
		
		// chr1:50|400_junction_exact
		this.expectedTranslations.add("SPYLAVSL");
		this.expectedTranslations.add("PLTLPSA");
		this.expectedTranslations.add("PLPCRQP");
		
		// chr1:200|400_junction_exact
		this.expectedTranslations.add("SPXLAVSL");
		this.expectedTranslations.add("PLNLPSA");
		this.expectedTranslations.add("PLTCRQP");
	}
}
