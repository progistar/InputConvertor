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

import zhanglab.inputconvertor.data.Exon;
import zhanglab.inputconvertor.data.GTFLoader;
import zhanglab.inputconvertor.data.GenomeLoader;
import zhanglab.inputconvertor.data.Transcript;
import zhanglab.inputconvertor.function.Translator;
import zhanglab.inputconvertor.input.CIRIquant;

public class RunTranslation {
	
	public static void test () throws IOException {
		File genomeFile = new File("/Users/seunghyukchoi/Documents/_resources/_databases/GRCh38.primary_assembly.genome.fa");
		File testFile = new File("/Users/seunghyukchoi/Documents/1_Projects/2023_Neoflow2/2_iRNAseq/CIRIquant/C3L-00973.T/CIRIquant_total_rnaseq.gtf");
		File referenceFile = new File("/Users/seunghyukchoi/Documents/_resources/_databases/gencode.v34.basic.annotation.gtf");
		
		GenomeLoader gmL = new GenomeLoader(genomeFile);
		GTFLoader gtfRef = new GTFLoader(referenceFile);
        
        BufferedWriter BW = new BufferedWriter(new FileWriter("/Users/seunghyukchoi/Documents/_resources/_databases/test.fa"));
        
        CIRIquant ciriQuant = new CIRIquant(testFile);
        ciriQuant.enrollGenomeSequence(gmL);
        ciriQuant.enrollReferenceGTF(gtfRef);
        ciriQuant.writeEntry(BW, null);
        BW.close();
	}
	
	public static void main(String[] args) throws IOException, ParseException {
		long startTime = System.currentTimeMillis();
		System.out.println("Translator v0.0.0");
		
		test();
		System.exit(1);
		Options options = new Options();
		
		// Options
		options.addOption("s", "stringtie", true, "GTF file path from StringTie");
		options.addOption("a", "arriba", true, "TSV file path from Arriba");
		options.addOption("c", "ciriquant", true, "GTF file path from CIRIQuant");
		options.addOption("i", "irfinder", true, "tsv file path from IRFinder");
		options.addOption("r", "gtf", true, "Reference GTF file");
		options.addOption("g", "fasta", true, "Genome fasta file");
		options.addOption("p", "pattern", true, "batch pattern");
		
		CommandLineParser parser = new DefaultParser();
        CommandLine cmd = parser.parse(options, args);
        
        
        
        // parse input
        //String mode = cmd.getOptionValue("m");
        String stringTieFilePath = cmd.getOptionValue("s");
        String genomeFastaFile = cmd.getOptionValue("g");
		
        File genomeFile = new File(genomeFastaFile);
        File stringTieFile = new File(stringTieFilePath);
        
        GenomeLoader gmL = new GenomeLoader(genomeFile);
        
        
        GTFLoader gtfL = new GTFLoader(stringTieFile);
        
        BufferedWriter BW = new BufferedWriter(new FileWriter("/Users/seunghyukchoi/Documents/_resources/_databases/test.fa"));
        
        gtfL.geneToTranscripts.forEach((g, ts)->{
        	StringBuilder sequence = new StringBuilder();
        	for(Transcript t : ts) {
        		sequence.setLength(0);
        		String protein = null;
        		for(Exon e : t.exons) {
        			sequence.append(gmL.getSequence(t.chr, e.start-1, e.end));
        		}
        		
        		try {
        			
        			for(int i=0; i<3; i++) {
        				if(t.strand.equalsIgnoreCase("+")) {
                			protein = Translator.translation(sequence.toString(), i);
                		} else {
                			protein = Translator.reverseComplementTranslation(sequence.toString(), i);
                		}
        				BW.append(">").append(t.tID).append("|").append("Frame_"+i);
            			BW.newLine();
            			BW.append(protein);
            			BW.newLine();
        			}
        		}catch(IOException ioe) {
        			
        		}
        	}
        });
        
        BW.close();
        
        long endTime = System.currentTimeMillis();
        
        System.out.println((endTime-startTime)/1000+" sec");
	}
}
