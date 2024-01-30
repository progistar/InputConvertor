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
import zhanglab.inputconvertor.data.VEPLoader;
import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.input.Arriba;
import zhanglab.inputconvertor.input.CIRIquant;
import zhanglab.inputconvertor.input.IRFinder;
import zhanglab.inputconvertor.input.StringTie;

public class RunTranslation {
 	
	public static ArrayList<String> newArgs = new ArrayList<String>();
	
	public static boolean checkArg (String short_, String long_, String arg1, String arg2) {
		boolean is = false;
		if((arg1.equalsIgnoreCase(short_) || arg1.equalsIgnoreCase(long_)) 
				&& arg2.length() != 0 && !arg2.startsWith("-")) {
			is = true;
			newArgs.add(arg1);
			newArgs.add(arg2);
		}
		return is;
	}
	
	public static void main(String[] args) throws IOException, ParseException {
		long startTime = System.currentTimeMillis();
		System.out.println(InputConvertorConstants.VERSION);
		Options options = new Options();
		
		// Options
		for(int i=0; i<args.length; i++) {
			if(i == args.length-1) {
				break;
			}
			
			
			if(checkArg("-g", "--genome",args[i], args[i+1])) {
				options.addOption("g", "genome", true, "Genome fasta file");
				System.out.println(i+"\t"+args[i]+"\t"+args[i+1]);
				
				
			} else if(checkArg("-p", "--protein",args[i], args[i+1])) {
				options.addOption("p", "protein", true, "Reference protein fasta file");
				System.out.println(i+"\t"+args[i]+"\t"+args[i+1]);
				
			} else if(checkArg("-r", "--gtf",args[i], args[i+1])) {
				options.addOption("r", "gtf", true, "Reference GTF file");
				System.out.println(i+"\t"+args[i]+"\t"+args[i+1]);
				
			} else if(checkArg("-o", "--output",args[i], args[i+1])) {
				options.addOption("o", "output", true, "Output path");
				System.out.println(i+"\t"+args[i]+"\t"+args[i+1]);
				
			} else if(checkArg("-s", "--stringtie",args[i], args[i+1])) {
				options.addOption("s", "stringtie", true, "GTF file path from StringTie");
				System.out.println(i+"\t"+args[i]+"\t"+args[i+1]);
				
			}  else if(checkArg("-f", "--fpkm",args[i], args[i+1])) {
				options.addOption("f", "fpkm", true, "FPKM threshold for StringTie (filter FPKM < value). Default 1.00");
				System.out.println(i+"\t"+args[i]+"\t"+args[i+1]);
				
			}  else if(checkArg("-a", "--arriba",args[i], args[i+1])) {
				options.addOption("a", "arriba", true, "TSV file path from Arriba");
				System.out.println(i+"\t"+args[i]+"\t"+args[i+1]);
				
			}  else if(checkArg("-c", "--ciriquant",args[i], args[i+1])) {
				options.addOption("c", "ciriquant", true, "GTF file path from CIRIQuant");
				System.out.println(i+"\t"+args[i]+"\t"+args[i+1]);
				
			}  else if(checkArg("-i", "--irfinder",args[i], args[i+1])) {
				options.addOption("i", "irfinder", true, "tsv file path from IRFinder");
				System.out.println(i+"\t"+args[i]+"\t"+args[i+1]);
				
			}  else if(checkArg("-v", "--vep_germ",args[i], args[i+1])) {
				options.addOption("v", "vep_germ", true, "Germline VEP file path from Ensembl Variant Effeect Predictor");
				System.out.println(i+"\t"+args[i]+"\t"+args[i+1]);
				
			}  else if(checkArg("-V", "--vep_soma",args[i], args[i+1])) {
				options.addOption("V", "vep_soma", true, "Somatic VEP file path from Ensembl Variant Effeect Predictor");
				System.out.println(i+"\t"+args[i]+"\t"+args[i+1]);
				
			}
		}
		
		
		args = new String[newArgs.size()];
		for(int i=0; i<args.length; i++) {
			args[i] = newArgs.get(i);
		}
		
		
		CommandLineParser parser = new DefaultParser();
        CommandLine cmd = parser.parse(options, args);
        
        File stringTieFile = null; File arribaFile = null; File ciriquantFile= null; File irfinderFile = null;
        File refGTFFile = null; File refProteinFile = null; String outputPath = null;
        File refGenomeFile = null; File vepSomaFile = null; File vepGermFile = null;  
        double fpkmThreshold = 1.00;
        
        if(cmd.hasOption("s")) {
        	stringTieFile = new File(cmd.getOptionValue("s"));
        	if(cmd.getOptionValue("f") != null) {
        		fpkmThreshold = Double.parseDouble(cmd.getOptionValue("f"));
        	}
        	System.out.println("FPKM < "+fpkmThreshold+" will be discarded in the given StringTie result");
        }
        if(cmd.hasOption("a")) {
        	arribaFile = new File(cmd.getOptionValue("a"));
        }
        if(cmd.hasOption("c")) {
        	ciriquantFile = new File(cmd.getOptionValue("c"));
        }
        if(cmd.hasOption("i")) {
        	irfinderFile = new File(cmd.getOptionValue("i"));
        }
        if(cmd.hasOption("r")) {
        	refGTFFile = new File(cmd.getOptionValue("r"));
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
        if(cmd.hasOption("g")) {
        	refGenomeFile = new File(cmd.getOptionValue("g"));
        }
        if(cmd.hasOption("o")) {
        	outputPath = cmd.getOptionValue("o");
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
        if(outputPath != null) {
        	System.out.println("Output sequence database: "+outputPath);
        } else {
        	System.out.println("Missing output path...");
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
        
        BufferedWriter BW = new BufferedWriter(new FileWriter(outputPath));
        
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
