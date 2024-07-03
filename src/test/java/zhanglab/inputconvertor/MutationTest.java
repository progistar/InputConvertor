package zhanglab.inputconvertor;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import org.junit.jupiter.api.Test;

import junit.framework.Assert;
import zhanglab.inputconvertor.data.FastaEntry;
import zhanglab.inputconvertor.data.GTFLoader;
import zhanglab.inputconvertor.data.GenomeLoader;
import zhanglab.inputconvertor.data.VARLoader;
import zhanglab.inputconvertor.input.Reference;

public class MutationTest {

	public File genomeFile = null;
	public File referenceFile = null;
	public File vepFile = null;
	public Hashtable<String, String> expectedTranslations = null;
	
	@Test
	public void runTest () throws IOException {

		MutationTest test = new MutationTest();
		test.loadTestData();
		
		GTFLoader refGTF = new GTFLoader(test.referenceFile);
		GenomeLoader gLoader = new GenomeLoader(test.genomeFile);
		VARLoader vepLoader = null;
		if(test.vepFile != null) {
			vepLoader = new VARLoader(test.vepFile);
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
        
        
        for(int i=0; i<entries.size(); i++) {
        	FastaEntry entry = entries.get(i);
        	String key = entry.description +"_frame"+entry.frame;
        	if(test.expectedTranslations.get(key) != null) {
        		Assert.assertTrue("Translated sequence was not equal: "+i +" "+ test.expectedTranslations.get(key) + "=>" + entry.sequence
        				, test.expectedTranslations.get(key).equalsIgnoreCase(entry.sequence));
        	}
        	
        }
        
		
	}
	
	public void loadTestData () {
		this.genomeFile = new File("test/test.genome.fa");
		this.referenceFile = new File("test/comb_mut_ref.gtf");
		this.vepFile = new File("test/comb_mut_vep.tsv");
		this.expectedTranslations = new Hashtable<String, String>();
		
		this.expectedTranslations.put("@chr1:1-1[SNP]G>T@chr1:1-1[INS]CCCC@chr1:2-4@chr1:5-5[SNP]A>C@chr1:6-20_frame0", "SPYLAVSL");
		this.expectedTranslations.put("@chr1:1-1[SNP]G>T@chr1:1-1[INS]CCCC@chr1:2-4@chr1:5-5[SNP]A>C@chr1:6-20_frame1", "PLTLPSA");
		this.expectedTranslations.put("@chr1:1-1[SNP]G>T@chr1:1-1[INS]CCCC@chr1:2-4@chr1:5-5[SNP]A>C@chr1:6-20_frame2", "PLPCRQP");
		this.expectedTranslations.put("@chr1:1-1[SNP]G>T@chr1:1-1[INS]CCCC@chr1:2-4@chr1:5-5@chr1:6-20_frame0", "SPXLAVSL");
		this.expectedTranslations.put("@chr1:1-1[SNP]G>T@chr1:1-1[INS]CCCC@chr1:2-4@chr1:5-5@chr1:6-20_frame1", "PLNLPSA");
		this.expectedTranslations.put("@chr1:1-1[SNP]G>T@chr1:1-1[INS]CCCC@chr1:2-4@chr1:5-5@chr1:6-20_frame2", "PLTCRQP");
	}
}
