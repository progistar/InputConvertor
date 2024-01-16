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
import zhanglab.inputconvertor.data.GTFLoader;
import zhanglab.inputconvertor.data.GenomeLoader;
import zhanglab.inputconvertor.data.Mutation;
import zhanglab.inputconvertor.data.VEPLoader;
import zhanglab.inputconvertor.input.Arriba;
import zhanglab.inputconvertor.input.CIRIquant;
import zhanglab.inputconvertor.input.StringTie;

public class RunTranslation {
	
	public static void testVEP() throws IOException {
		File vepFile = new File("/Users/seunghyukchoi/Documents/1_Projects/2023_Neoflow2/2_iRNAseq/vep/C3L-01632.txt");
		VEPLoader vep = new VEPLoader(vepFile, true);
		
		ArrayList<Mutation> mutations = vep.getDELByRange("chr19", 9126193, 9126195);
		
		for(Mutation mutation : mutations) {
			System.out.println(mutation.toString());
		}
		
		
	}
	
 	public static void testAll () throws IOException {
		File genomeFile = new File("/Users/seunghyukchoi/Documents/_resources/_databases/GRCh38.primary_assembly.genome.fa");
		File testFile = new File("/Users/seunghyukchoi/Documents/1_Projects/2023_Neoflow2/2_iRNAseq/CIRIquant/C3L-00973.T/CIRIquant_total_rnaseq.gtf");
		File referenceFile = new File("/Users/seunghyukchoi/Documents/_resources/_databases/gencode.v34.basic.annotation.gtf");
		File arribaFile = new File("/Users/seunghyukchoi/Documents/1_Projects/2023_Neoflow2/2_iRNAseq/arriba/C3L-00973.T.arriba.fusions.tsv");
		
		GenomeLoader gmL = new GenomeLoader(genomeFile);
		GTFLoader gtfRef = new GTFLoader(referenceFile);
        
        
        CIRIquant ciriQuant = new CIRIquant(testFile);
        ciriQuant.enrollGenomeSequence(gmL);
        ciriQuant.enrollReferenceGTF(gtfRef);
        ArrayList<FastaEntry> entries = ciriQuant.getFastaEntry();
        
        Arriba arriba = new Arriba(arribaFile);
        entries.addAll(arriba.getFastaEntry());
        
        BufferedWriter BW = new BufferedWriter(new FileWriter("/Users/seunghyukchoi/Documents/_resources/_databases/test.fa"));
        
        for(FastaEntry entry : entries) {
        	BW.append(">"+entry.tool+"|"+entry.header);
        	BW.newLine();
        	BW.append(entry.sequence);
        	BW.newLine();
        }
        
        BW.close();
	}
 	
 	public static void testStringTie () throws IOException {
 		File genomeFile = new File("/Users/seunghyukchoi/Documents/_resources/_databases/GRCh38.primary_assembly.genome.fa");
		File stringTieFile = new File("/Users/seunghyukchoi/Documents/1_Projects/2023_Neoflow2/2_iRNAseq/stringtie/C3L-00973.T.stringtie_output.gtf");
		GenomeLoader gmL = new GenomeLoader(genomeFile);
		StringTie stringTie = new StringTie(stringTieFile);
		stringTie.enrollGenomeSequence(gmL);
		
		ArrayList<FastaEntry> entries = stringTie.getFastaEntry();
		
		BufferedWriter BW = new BufferedWriter(new FileWriter("/Users/seunghyukchoi/Documents/_resources/_databases/test.fa"));
        
        for(FastaEntry entry : entries) {
        	BW.append(">"+entry.tool+"|"+entry.header);
        	BW.newLine();
        	BW.append(entry.sequence);
        	BW.newLine();
        }
        
        BW.close();
 	}
	
	public static void main(String[] args) throws IOException, ParseException {
		long startTime = System.currentTimeMillis();
		System.out.println("Translator v0.0.0");
		
		testStringTie();
		System.exit(1);
		Options options = new Options();
		
		// Options
		options.addOption("s", "stringtie", true, "GTF file path from StringTie");
		options.addOption("a", "arriba", true, "TSV file path from Arriba");
		options.addOption("c", "ciriquant", true, "GTF file path from CIRIQuant");
		options.addOption("i", "irfinder", true, "tsv file path from IRFinder");
		options.addOption("r", "gtf", true, "Reference GTF file");
		options.addOption("v", "VEP", true, "VEP file path from Ensembl Variant Effeect Predictor");
		options.addOption("g", "fasta", true, "Genome fasta file");
		options.addOption("p", "pattern", true, "batch pattern");
		
		CommandLineParser parser = new DefaultParser();
        CommandLine cmd = parser.parse(options, args);
        
        
        
        // parse input
        //String mode = cmd.getOptionValue("m");
        String stringTieFilePath = cmd.getOptionValue("s");
        String genomeFastaFile = cmd.getOptionValue("g");
		
        long endTime = System.currentTimeMillis();
        
        System.out.println((endTime-startTime)/1000+" sec");
	}
}
