package zhanglab.inputconvertor.run;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.data.FastaEntry;
import zhanglab.inputconvertor.data.FastaLoader;
import zhanglab.inputconvertor.data.GTFLoader;
import zhanglab.inputconvertor.data.GenomeLoader;
import zhanglab.inputconvertor.data.VEPLoader;
import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.input.Arriba;
import zhanglab.inputconvertor.input.CIRIquant;
import zhanglab.inputconvertor.input.IRFinder;
import zhanglab.inputconvertor.input.Reference;
import zhanglab.inputconvertor.input.StringTie;

public class RunTranslation {
 	
	public static File stringTieFile = null; 
	public static File arribaFile = null; 
	public static File ciriquantFile= null; 
	public static File irfinderFile = null;
	public static File refGTFFile = null; 
	public static File refProteinFile = null; 
	public static File outputFile = null;
	public static File tableFile = null;
	public static File refGenomeFile = null; 
	public static File vepSomaFile = null; 
	public static File vepGermFile = null;
	public static String uniqueId = null;
    public static double fpkmThreshold = 1.00;
	
	
	public static void main(String[] args) throws IOException, ParseException {
		parseOptions(args);
		
		long startTime = System.currentTimeMillis();
		System.out.println(InputConvertorConstants.VERSION);
		
        
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
        Reference reference = new Reference(refGTF);
        reference.enrollGenomeSequence(gmL);
        entries.addAll(reference.getFastaEntry());
        
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
        
        BufferedWriter BW = new BufferedWriter(new FileWriter(outputFile));
        BufferedWriter table = new BufferedWriter(new FileWriter(tableFile));
        
        table.append("EntryId\tGeneId\tTranscriptId\tFrame\tStrand\tExons");
        table.newLine();
        for(FastaEntry entry : entries) {
        	BW.append(">"+entry.toHeader(uniqueId));
        	BW.newLine();
        	BW.append(entry.sequence);
        	BW.newLine();
        	
        	table.append(entry.toMeta(uniqueId));
        	table.newLine();
        }
        
        BW.close();
        table.close();
        
        long endTime = System.currentTimeMillis();
        
        System.out.println((endTime-startTime)/1000+" sec");
	}
	
	
	public static void parseOptions (String[] args) {
		CommandLine cmd = null;
		Options options = new Options();
		
		// Mandatory
		Option optionGenome = Option.builder("g")
				.longOpt("genome").argName("fasta")
				.hasArg()
				.required(true)
				.desc("reference genome sequence file")
				.build();
		
		Option optionProtein = Option.builder("p")
				.longOpt("protein").argName("fasta")
				.hasArg()
				.required(false)
				.desc("reference protein sequence file")
				.build();
		
		Option optionGTF = Option.builder("r")
				.longOpt("gene_model").argName("gtf")
				.hasArg()
				.required(true)
				.desc("reference gtf file")
				.build();
		
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("fasta")
				.hasArg()
				.required(true)
				.desc("output file")
				.build();
		
		Option optionTable = Option.builder("t")
				.longOpt("table").argName("tsv")
				.hasArg()
				.required(true)
				.desc("output table file")
				.build();
		
		Option optionUniqueId = Option.builder("u")
				.longOpt("uid").argName("string")
				.hasArg()
				.required(true)
				.desc("unique id for the file")
				.build();
		
		Option optionStringTie = Option.builder("s")
				.longOpt("stringtie").argName("gtf")
				.hasArg()
				.required(false)
				.desc("StringTie GTF file")
				.build();
		
		Option optionFPKM = Option.builder("f")
				.longOpt("fpkm").argName("float")
				.hasArg()
				.required(false)
				.desc("FPKM cutoff value in StringTie")
				.build();
		
		Option optionArriba = Option.builder("a")
				.longOpt("arriba").argName("tsv")
				.hasArg()
				.required(false)
				.desc("Arriba TSV file")
				.build();
		
		Option optionCIRIquant = Option.builder("c")
				.longOpt("ciriquant").argName("gtf")
				.hasArg()
				.required(false)
				.desc("CIRIquant GTF file")
				.build();
		
		Option optionIRFinder = Option.builder("i")
				.longOpt("irfinder").argName("tsv")
				.hasArg()
				.required(false)
				.desc("IRFinder TSV file")
				.build();
		
		Option optionVEPGerm = Option.builder("v")
				.longOpt("germ").argName("tsv")
				.hasArg()
				.required(false)
				.desc("VEP TSV file")
				.build();
		
		Option optionVEPSoma = Option.builder("V")
				.longOpt("soma").argName("tsv")
				.hasArg()
				.required(false)
				.desc("VEP TSV file")
				.build();
		
		
		options.addOption(optionGenome)
		.addOption(optionProtein)
		.addOption(optionGTF)
		.addOption(optionStringTie)
		.addOption(optionFPKM)
		.addOption(optionArriba)
		.addOption(optionCIRIquant)
		.addOption(optionIRFinder)
		.addOption(optionVEPGerm)
		.addOption(optionVEPSoma)
		.addOption(optionOutput)
		.addOption(optionTable)
		.addOption(optionUniqueId);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
	    ArrayList<String> tmpArgs = new ArrayList<String>();
	    for(int i=0; i<args.length; i++) {
	    	if(args[i].equalsIgnoreCase("-g") || args[i].equalsIgnoreCase("--genome") ||
			args[i].equalsIgnoreCase("-p") || args[i].equalsIgnoreCase("--protein") ||
			args[i].equalsIgnoreCase("-r") || args[i].equalsIgnoreCase("--gene_model") ||
			args[i].equalsIgnoreCase("-s") || args[i].equalsIgnoreCase("--stringtie") ||
			args[i].equalsIgnoreCase("-f") || args[i].equalsIgnoreCase("--fpkm") ||
			args[i].equalsIgnoreCase("-a") || args[i].equalsIgnoreCase("--arriba") ||
			args[i].equalsIgnoreCase("-c") || args[i].equalsIgnoreCase("--ciriquant") ||
			args[i].equalsIgnoreCase("-i") || args[i].equalsIgnoreCase("--irfinder") ||
			args[i].equalsIgnoreCase("-v") || args[i].equalsIgnoreCase("--germ") ||
			args[i].equalsIgnoreCase("-V") || args[i].equalsIgnoreCase("--soma") ||
			args[i].equalsIgnoreCase("-o") || args[i].equalsIgnoreCase("--output") ||
			args[i].equalsIgnoreCase("-t") || args[i].equalsIgnoreCase("--table") ||
			args[i].equalsIgnoreCase("-u") || args[i].equalsIgnoreCase("--uid")) {
	    		tmpArgs.add(args[i++]);
	    		tmpArgs.add(args[i]);
	    	}
	    }
	    
	    String[] nArgs = new String[tmpArgs.size()];
	    for(int i =0; i<tmpArgs.size(); i++) {
	    	nArgs[i] = tmpArgs.get(i);
	    }
	    
	    
		try {
		    cmd = parser.parse(options, nArgs, true);
		    
		    if(cmd.hasOption("g")) {
		    	refGenomeFile = new File(cmd.getOptionValue("g"));
		    }
		    
		    if(cmd.hasOption("p")) {
		    	refProteinFile = new File(cmd.getOptionValue("p"));
		    }
		    
		    if(cmd.hasOption("r")) {
		    	refGTFFile = new File(cmd.getOptionValue("r"));
		    }
		    
		    if(cmd.hasOption("s")) {
		    	stringTieFile = new File(cmd.getOptionValue("s"));
		    }
		    
		    if(cmd.hasOption("f")) {
		    	fpkmThreshold = Double.parseDouble(cmd.getOptionValue("f"));
		    }
		    
		    if(cmd.hasOption("a")) {
		    	arribaFile = new File(cmd.getOptionValue("a"));
		    }
		    
		    if(cmd.hasOption("c")) {
		    	ciriquantFile = new File(cmd.getOptionValue("c"));
		    }
		    
		    if(cmd.hasOption("v")) {
		    	vepGermFile = new File(cmd.getOptionValue("v"));
		    }
		    
		    if(cmd.hasOption("V")) {
		    	vepSomaFile = new File(cmd.getOptionValue("V"));
		    }
		    
		    if(cmd.hasOption("p")) {
		    	refProteinFile = new File(cmd.getOptionValue("p"));
		    }
		    
		    if(cmd.hasOption("o")) {
		    	outputFile = new File(cmd.getOptionValue("o"));
		    }
		    
		    if(cmd.hasOption("t")) {
		    	tableFile = new File(cmd.getOptionValue("t"));
		    }
		    
		    if(cmd.hasOption("u")) {
		    	uniqueId = cmd.getOptionValue("u");
		    }
		    
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			isFail = true;
		}
		
		if(isFail) {
		    helper.printHelp("Usage:", options);
		    System.exit(0);
		}
		
		System.out.println();
	}
}
