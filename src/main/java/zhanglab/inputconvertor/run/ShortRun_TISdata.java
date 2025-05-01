package zhanglab.inputconvertor.run;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.data.FastaEntry;
import zhanglab.inputconvertor.data.FastaLoader;
import zhanglab.inputconvertor.data.GTFLoader;
import zhanglab.inputconvertor.data.GenomeLoader;
import zhanglab.inputconvertor.data.VARLoader;
import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.input.Arriba;
import zhanglab.inputconvertor.input.Reference;
import zhanglab.inputconvertor.input.StringTie;

public class ShortRun_TISdata {

	public static File stringTieFile = null; 
	public static File arribaFile = null; 
	public static File ciriquantFile= null; 
	public static File irfinderFile = null;
	public static File refGTFFile = null; 
	public static File refProteinFile = null; 
	public static File outputFile = null;
	public static File refGenomeFile = null; 
	public static File varFile = null; 
	public static String uniqueId = null;
    public static double fpkmThreshold = 0.00;
	
	
	public static void main(String[] args) throws IOException, ParseException {
		
		long startTime = System.currentTimeMillis();
		System.out.println(InputConvertorConstants.VERSION);
		
        
        // requirements
		// stringTieFile = new File("/Volumes/Papers/2025_Emory/1.StringTieMerge/K21_D15_19_F1_final.annotated.gtf");
		refGenomeFile = new File("/Users/seunghyukchoi/Documents/_resources/_databases/GRCh38.primary_assembly.genome.fa");
		refGTFFile = new File("/Users/seunghyukchoi/Documents/_resources/_databases/gencode.v42.basic.annotation.gtf");
		outputFile = new File("/Volumes/Papers/2025_TIS/Initial_data.dfasta");
		uniqueId = "";
		
        if(stringTieFile != null) {
        	System.out.println("Novel isoforms from StringTie: " + stringTieFile.getName());
        }
        if(arribaFile != null) {
        	System.out.println("Fusion genes from Arriba: " + arribaFile.getName());
        }
        if(ciriquantFile != null) {
        	System.out.println("CircRNAs from CIRIquant: " + ciriquantFile.getName());
        }
        if(irfinderFile != null) {
        	System.out.println("Reteined-introns from IRFinder: " + irfinderFile.getName());
        }
        if(varFile != null) {
        	System.out.println("Variant calls from: " + varFile.getName());
        }
        if(refGTFFile != null) {
        	System.out.println("Reference transcriptome model: " + refGTFFile.getName());
        }
        if(refGenomeFile != null) {
        	System.out.println("Reference genome sequences: " + refGenomeFile.getName());
        }
        if(refProteinFile != null) {
        	System.out.println("Reference protein database: " + refProteinFile.getName());
        	System.out.println("Reference protein sequences will be appended to at the end of a customized database.");
        }
        
        
        // load GML
        GenomeLoader gmL = new GenomeLoader(refGenomeFile);
        if(varFile != null) {
        	VARLoader var = new VARLoader(varFile);
        	gmL.enrollVEPLaoder(var);
        }
        GTFLoader refGTF = new GTFLoader(refGTFFile);
        
        // Do StringTie
        ArrayList<FastaEntry> entries = new ArrayList<FastaEntry>();

        // Do reference fasta
        Reference reference = new Reference(refGTF);
        reference.enrollGenomeSequence(gmL);
        //entries.addAll(reference.getFastaEntry(true));
        entries.addAll(reference.getFastaEntry(false));
        
        if(refProteinFile != null) {
        	FastaLoader refProteins = new FastaLoader(refProteinFile);
        	entries.addAll(refProteins.entries);
        }
        
        if(stringTieFile != null) {
        	StringTie stringTie = new StringTie(stringTieFile);
        	stringTie.enrollGenomeSequence(gmL);
        	entries.addAll(stringTie.getFastaEntry(fpkmThreshold));
        }
        
        // Do Arriba
        if(arribaFile != null) {
        	Arriba arriba = new Arriba(arribaFile);
            entries.addAll(arriba.getFastaEntry());
        }
        
        
        System.out.println("A total of entries: "+entries.size());
        
        BufferedWriter BW = new BufferedWriter(new FileWriter(outputFile));
        
        // EntryId|GeneId|TranscriptId|GeneName|Frame|Strand|Exons
        for(FastaEntry entry : entries) {
        	BW.append(">"+entry.toMeta(""));
    		BW.newLine();
    		BW.append(entry.sequence);
    		BW.newLine();
    		BW.append(entry.nucleotide);
    		BW.newLine();
        }
        
        BW.close();
        
        long endTime = System.currentTimeMillis();
        
        System.out.println((endTime-startTime)/1000+" sec");
	}
}
