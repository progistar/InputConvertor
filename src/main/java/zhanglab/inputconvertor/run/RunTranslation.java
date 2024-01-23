package zhanglab.inputconvertor.run;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.data.FastaEntry;
import zhanglab.inputconvertor.data.FastaLoader;
import zhanglab.inputconvertor.data.GTFLoader;
import zhanglab.inputconvertor.data.GenomeLoader;
import zhanglab.inputconvertor.data.Mutation;
import zhanglab.inputconvertor.data.VEPLoader;
import zhanglab.inputconvertor.input.Arriba;
import zhanglab.inputconvertor.input.CIRIquant;
import zhanglab.inputconvertor.input.IRFinder;
import zhanglab.inputconvertor.input.StringTie;

public class RunTranslation {
 	
	public static void main(String[] args) throws IOException, ParseException {
		long startTime = System.currentTimeMillis();
		System.out.println("Translator v0.0.1");
		Options options = new Options();
		
		// Options
		options.addOption("s", "stringtie", true, "GTF file path from StringTie");
		options.addOption("f", "fpkm", true, "FPKM threshold for StringTie (filter FPKM < value). Default 1.00");
		options.addOption("a", "arriba", true, "TSV file path from Arriba");
		options.addOption("c", "ciriquant", true, "GTF file path from CIRIQuant");
		options.addOption("i", "irfinder", true, "tsv file path from IRFinder");
		options.addOption("r", "gtf", true, "Reference GTF file");
		options.addOption("v", "VEP_germ", true, "Germline VEP file path from Ensembl Variant Effeect Predictor");
		options.addOption("V", "VEP_soma", true, "Somatic VEP file path from Ensembl Variant Effeect Predictor");
		options.addOption("g", "fasta", true, "Genome fasta file");
		options.addOption("p", "fasta", true, "Reference protein fasta file");
		options.addOption("o", "output_prefix", true, "output_prefix");
		
		CommandLineParser parser = new DefaultParser();
        CommandLine cmd = parser.parse(options, args);
        
        File stringTieFile = null; File arribaFile = null; File ciriquantFile= null; File irfinderFile = null;
        File refGTFFile = null; File refProteinFile = null; String outputPrefix = null;
        File refGenomeFile = null; File vepSomaFile = null; File vepGermFile = null;  
        double fpkmThreshold = 1.00;
        
        if(cmd.getOptionValue("s") != null) {
        	stringTieFile = new File(cmd.getOptionValue("s"));
        	if(cmd.getOptionValue("f") != null) {
        		fpkmThreshold = Double.parseDouble(cmd.getOptionValue("f"));
        	}
        	System.out.println("FPKM < "+fpkmThreshold+" will be discarded in the given StringTie result");
        }
        if(cmd.getOptionValue("a") != null) {
        	arribaFile = new File(cmd.getOptionValue("a"));
        }
        if(cmd.getOptionValue("c") != null) {
        	ciriquantFile = new File(cmd.getOptionValue("c"));
        }
        if(cmd.getOptionValue("i") != null) {
        	irfinderFile = new File(cmd.getOptionValue("i"));
        }
        if(cmd.getOptionValue("r") != null) {
        	refGTFFile = new File(cmd.getOptionValue("r"));
        }
        if(cmd.getOptionValue("v") != null) {
        	vepGermFile = new File(cmd.getOptionValue("v"));
        }
        if(cmd.getOptionValue("V") != null) {
        	vepSomaFile = new File(cmd.getOptionValue("V"));
        }
        if(cmd.getOptionValue("p") != null) {
        	refProteinFile = new File(cmd.getOptionValue("p"));
        }
        if(cmd.getOptionValue("g") != null) {
        	refGenomeFile = new File(cmd.getOptionValue("g"));
        }
        if(cmd.getOptionValue("o") != null) {
        	outputPrefix = cmd.getOptionValue("o");
        }
        
        // requirements
        
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
        if(vepGermFile != null) {
        	System.out.println("Germline variant calls from: " + vepGermFile.getName());
        }
        if(vepSomaFile != null) {
        	System.out.println("Somatic variant calls from: " + vepSomaFile.getName());
        }
        if(refGTFFile != null) {
        	System.out.println("Reference transcriptome model: " + refGTFFile.getName());
        } else {
        	System.out.println("Missing reference transcriptome model...");
        	System.exit(1);
        }
        if(refGenomeFile != null) {
        	System.out.println("Reference genome sequences: " + refGenomeFile.getName());
        } else {
        	System.out.println("Missing reference genome sequences...");
        	System.exit(1);
        }
        if(refProteinFile != null) {
        	System.out.println("Reference protein database: " + refProteinFile.getName());
        	System.out.println("Reference protein sequences will be appended to at the end of a customized database.");
        }
        
        
        // load GML
        GenomeLoader gmL = new GenomeLoader(refGenomeFile);
        GTFLoader refGTF = new GTFLoader(refGTFFile);
        if(vepGermFile != null) {
        	VEPLoader vep = new VEPLoader(vepGermFile, false);
        	gmL.enrollGermVEPLaoder(vep);
        }
        if(vepSomaFile != null) {
        	VEPLoader vep = new VEPLoader(vepSomaFile, true);
        	gmL.enrollSomaticVEPLaoder(vep);
        }
        
        // Do StringTie
        ArrayList<FastaEntry> entries = new ArrayList<FastaEntry>();

        // Do reference fasta
        if(refProteinFile != null) {
        	FastaLoader refProteins = new FastaLoader(refProteinFile);
        	entries.addAll(refProteins.entries);
        }
        
        if(stringTieFile != null) {
        	StringTie stringTie = new StringTie(stringTieFile);
        	stringTie.enrollGenomeSequence(gmL);
        	entries.addAll(stringTie.getFastaEntry(fpkmThreshold));
        }
        
        // Do CIRIquant
        if(ciriquantFile != null) {
        	CIRIquant ciriquant = new CIRIquant(ciriquantFile);
        	ciriquant.enrollGenomeSequence(gmL);
        	ciriquant.enrollReferenceGTF(refGTF);
        	entries.addAll(ciriquant.getFastaEntry());
        }
        
        // Do IRFinder
        if(irfinderFile != null) {
        	IRFinder irfinder = new IRFinder(irfinderFile);
        	irfinder.enrollGenomeSequence(gmL);
        	irfinder.enrollReferenceGTF(refGTF);
        	entries.addAll(irfinder.getFastaEntry());
        }
        
        // Do Arriba
        if(arribaFile != null) {
        	Arriba arriba = new Arriba(arribaFile);
            entries.addAll(arriba.getFastaEntry());
        }
        
        
        System.out.println("A total of entries: "+entries.size());
        
        String outputFile = outputPrefix +".combined.fasta";
        BufferedWriter BW = new BufferedWriter(new FileWriter(outputFile));
        
        for(FastaEntry entry : entries) {
        	BW.append(">"+entry.toHeader());
        	BW.newLine();
        	BW.append(entry.sequence);
        	BW.newLine();
        }
        
        BW.close();
        
        
        long endTime = System.currentTimeMillis();
        
        System.out.println((endTime-startTime)/1000+" sec");
	}
}
