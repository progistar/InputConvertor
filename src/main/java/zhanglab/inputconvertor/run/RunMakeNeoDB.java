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
import zhanglab.inputconvertor.data.VARLoader;
import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.input.Arriba;
import zhanglab.inputconvertor.input.MANTA;
import zhanglab.inputconvertor.input.Reference;
import zhanglab.inputconvertor.input.StringTie;

public class RunMakeNeoDB {
 	
	public static File stringTieFile = null; 
	public static File arribaFile = null; 
	public static File refGTFFile = null; 
	public static File refProteinFile = null; 
	public static File outputFile = null;
	public static File refGenomeFile = null; 
	public static File varFile = null; 
	public static File MANTAFile = null;
	public static String uniqueId = null;
    public static double fpkmThreshold = 1.00;
    
	// ## Translation of mutated sequences
	public static int maxFlankLength	=	15;
	
	public static void main(String[] args) throws IOException, ParseException {
		System.out.println(InputConvertorConstants.VERSION);
		parseOptions(args);
		
		long startTime = System.currentTimeMillis();
		
        
        // requirements
        
        if(stringTieFile != null) {
        	System.out.println("Novel isoforms from StringTie: " + stringTieFile.getName());
        }
        if(arribaFile != null) {
        	System.out.println("Fusion genes from Arriba: " + arribaFile.getName());
        }
        if(varFile != null) {
        	System.out.println("Variant calls from: " + varFile.getName());
        }
        if(MANTAFile != null) {
        	System.out.println("Structural variations from : " + MANTAFile.getName());
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
        
        // Do StringTie
        ArrayList<FastaEntry> entries = new ArrayList<FastaEntry>();
        
        
        if(refGTFFile != null) {
        	GTFLoader refGTF = new GTFLoader(refGTFFile);
        	Reference reference = new Reference(refGTF);
        	reference.enrollGenomeSequence(gmL);
        	//entries.addAll(reference.getFastaEntry(true));
        	entries.addAll(reference.getFastaEntry(false));
        }
        
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
        
        // Do MANTA
        if(MANTAFile != null) {
        	MANTA manta = new MANTA(MANTAFile);
            entries.addAll(manta.getFastaEntry());
        }
        
        
        System.out.println("A total of entries: "+entries.size());
        
        BufferedWriter BW = new BufferedWriter(new FileWriter(outputFile));
        
        // EntryId|GeneId|TranscriptId|GeneName|Frame|Strand|Exons
        for(FastaEntry entry : entries) {
        	BW.append(">"+entry.toMeta(uniqueId));
        	BW.newLine();
        	BW.append(entry.sequence);
        	BW.newLine();
        }
        
        BW.close();
        
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
				.required(false)
				.desc("reference gtf file")
				.build();
		
		Option optionOutput= Option.builder("o")
				.longOpt("output").argName("string")
				.hasArg()
				.required(true)
				.desc("prefix of output file")
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
		
		Option optionVEP = Option.builder("v")
				.longOpt("var").argName("tsv|vcf")
				.hasArg()
				.required(false)
				.desc("VEP or VCF file")
				.build();
		
		Option optionMANTA = Option.builder("m")
				.longOpt("manta").argName("vcf")
				.hasArg()
				.required(false)
				.desc("VCF file from MANTA")
				.build();
		
		
		options.addOption(optionGenome)
		.addOption(optionProtein)
		.addOption(optionGTF)
		.addOption(optionStringTie)
		.addOption(optionFPKM)
		.addOption(optionArriba)
		.addOption(optionVEP)
		.addOption(optionOutput)
		.addOption(optionMANTA);
		
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
			args[i].equalsIgnoreCase("-v") || args[i].equalsIgnoreCase("--var") ||
			args[i].equalsIgnoreCase("-m") || args[i].equalsIgnoreCase("--manta") ||
			args[i].equalsIgnoreCase("-o") || args[i].equalsIgnoreCase("--output")) {
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
		    
		    if(cmd.hasOption("v")) {
		    	varFile = new File(cmd.getOptionValue("v"));
		    }
		    
		    if(cmd.hasOption("p")) {
		    	refProteinFile = new File(cmd.getOptionValue("p"));
		    }
		    
		    if(cmd.hasOption("m")) {
		    	MANTAFile = new File(cmd.getOptionValue("m"));
		    }
		    
		    if(cmd.hasOption("o")) {
		    	outputFile = new File(cmd.getOptionValue("o"));
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
