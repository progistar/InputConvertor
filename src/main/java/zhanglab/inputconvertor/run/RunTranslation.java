package zhanglab.inputconvertor.run;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

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

public class RunTranslation {
	
	public static void main(String[] args) throws IOException, ParseException {
		long startTime = System.currentTimeMillis();
		System.out.println("Translator v0.0.0");
		
		Options options = new Options();
		
		// Options
		options.addOption("i", "in", true, "input file");
		options.addOption("t", "tool", true, "input tool (ex> stringtie, ciriquant, arriba, irfinder)");
		options.addOption("g", "gtf", true, "reference GTF file");
		options.addOption("f", "fasta", true, "genome fasta file");
		
		CommandLineParser parser = new DefaultParser();
        CommandLine cmd = parser.parse(options, args);
        
        // parse input
        //String mode = cmd.getOptionValue("m");
        //String inputTool = cmd.getOptionValue("t");
        //String genomeFastaFile = cmd.getOptionValue("f");
		
        GenomeLoader gmL = new GenomeLoader(new File("/Users/seunghyukchoi/Documents/_resources/_databases/GRCh38.primary_assembly.genome.fa"));
        GTFLoader gtfL = new GTFLoader(new File("/Users/seunghyukchoi/Documents/_resources/_databases/gencode.v42.basic.annotation.gtf"));
        
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
