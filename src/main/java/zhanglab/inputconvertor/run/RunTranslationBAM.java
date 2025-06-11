package zhanglab.inputconvertor.run;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import zhanglab.inputconvertor.data.FastaEntry;
import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.function.Mode;
import zhanglab.inputconvertor.function.Translator;

public class RunTranslationBAM {

	public static File bamFile = null;
	public static File outputFile = null;
	public static int maxMer = 25;
	
	
	public static void main(String[] args) throws IOException {
		parseOptions(args);
		
		System.out.println("## Unmapped translation ##");
		System.out.println("Load "+bamFile.getName());
		writeFastaEntry();
	}
	
	public static void writeFastaEntry () throws IOException {
		ArrayList<FastaEntry> fastaEntries = new ArrayList<FastaEntry>();
		
		String strandedness = "NON";
		
		int R1F = 0;
		int R1R = 0;
		int R2F = 0;
		int R2R = 0;
		
		try (SamReader samReader = SamReaderFactory.makeDefault().open(bamFile)) {
			SAMRecordIterator iterator = samReader.iterator();
			int size = 1000000;
			while((size--) > 0 && iterator.hasNext()) {
				SAMRecord samRecord = iterator.next();
				
				boolean isPass = false;
				if(samRecord.isSecondaryAlignment()) {
            		isPass = true;
            	}
				
				Object xsTag = samRecord.getAttribute("XS");
				if(xsTag == null) {
					isPass = true;
				}
				
				if(isPass) {
					continue;
				}
				
				
				int flags = samRecord.getFlags();
				boolean isFirstSegment = (0x40 & flags) == 0x40 ? true : false;
				boolean isForward = (0x10 & flags) == 0x10 ? false : true;
				boolean strand = ((Character) xsTag) == '+' ? true : false;
				
				// first segment
				if(isFirstSegment) {
					if(isForward == strand) {
						R1F++;
					} else {
						R1R++;
					}
				} 
				// second segment
				else {
					if(isForward == strand) {
						R2F++;
					} else {
						R2R++;
					}
				}
				
			}
			iterator.close();
			
			if(R1F*10 < R1R && R2F > R2R*10) {
				strandedness = InputConvertorConstants.RF_STRANDED;
			} else if(R1F > R1R*10 && R2F*10 < R2R) {
				strandedness = InputConvertorConstants.FR_STRANDED;
			} 
			// for single end
			else if( (R1F + R2F) > 10 * (R1R + R2R) ) {
				strandedness = InputConvertorConstants.F_STRANDED;
			} else if( 10 * (R1F + R2F) < (R1R + R2R) ) {
				strandedness = InputConvertorConstants.R_STRANDED;
			}
			// not found
			else {
				strandedness = InputConvertorConstants.NON_STRANDED;
			}
			
			System.out.println("Estimate strandedness");
			System.out.println("1F\t1R\t2F\t2R");
			System.out.println(R1F+"\t"+R1R+"\t"+R2F+"\t"+R2R);
			
			if(R1F+R1R+R2F+R2R == 0) {
				System.out.println("Fail to estimate stradedness!");
				System.out.println("It looks single-end RNA-seq experiement. Please specify strandedness.");
				System.exit(1);
			} else {
				System.out.println(strandedness);
			}
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		try (SamReader samReader = SamReaderFactory.makeDefault().open(bamFile)) {
			SAMRecordIterator iterator = samReader.queryUnmapped();
			int totalUmap = 0;
			while (iterator.hasNext()) {
				iterator.next();
				totalUmap++;
			}
			System.out.println("Total unmapped reads: "+totalUmap);
			iterator.close();
		}
		
		try (SamReader samReader = SamReaderFactory.makeDefault().open(bamFile)) {
			SAMRecordIterator iterator = samReader.queryUnmapped();
			Hashtable<String, Integer> readMap = new Hashtable<String, Integer>();
			
			while (iterator.hasNext()) {
	            SAMRecord samRecord = iterator.next();
	            int flags = samRecord.getFlags();
	            ArrayList<Character> strands = Mode.getStrandedness(flags, strandedness);
	            
	            for(Character strand : strands) {
	            	String sequence = null;
	            	if(strand == '+') {
	            		sequence = samRecord.getReadString();
	            	} else {
	            		sequence = Translator.getReverseComplement(samRecord.getReadString());
	            	}
	            	
	            	for(int fr=0; fr<3; fr++) {
            			String peptide = Translator.translation(sequence, fr)[1];
            			int len = peptide.length();
            			for(int i=0; i<len; i++) {
            				int start = i;
            				int end = start + Math.min(maxMer, len);
            				
            				if(end > len) {
            					break;
            				} else {
            					String subPeptide = peptide.substring(start, end);
            					String[] xSplit = subPeptide.split("X");
            					
            					for(String p : xSplit) {
            						if(p.length() < InputConvertorConstants.MIN_PEPT_LEN) {
            							continue;
            						}
            						
            						Integer reads = readMap.get(p);
            						if(reads == null) {
            							reads = 0;
            						}
            						
            						readMap.put(p, reads+1);
            					}
            				}
            			}
            		}
	            }
			}
			
			iterator.close();
			
			int[] readCounts = new int[101];
			readMap.forEach((pept, reads)->{
				int idx = Math.min(reads, readCounts.length-1);
				readCounts[idx]++;
				if(reads >= 10) {
					FastaEntry entry = new FastaEntry();
					entry.tool = InputConvertorConstants.UNMAPPED_HEADER_ID;
					entry.idx = fastaEntries.size()+1;
					entry.sequence = pept;
					entry.strand = ".";
					entry.frame = ".";
					entry.transcriptId = ".";
					entry.geneId = ".";
					entry.geneName = ".";
					entry.description = "Unmapped";
					fastaEntries.add(entry);
				}
			});
			
			System.out.println("Length\tCount");
			for(int i=1; i<readCounts.length; i++) {
				System.out.println(i+"\t"+readCounts[i]);
			}
			
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		
		 BufferedWriter BW = new BufferedWriter(new FileWriter(outputFile));
	        
	        // EntryId|GeneId|TranscriptId|GeneName|Frame|Strand|Exons
	        for(FastaEntry entry : fastaEntries) {
	        	BW.append(">"+entry.toMeta(null));
	        	BW.newLine();
	        	BW.append(entry.sequence);
	        	BW.newLine();
	        }
	        
	        BW.close();
	}
	
	

	public static void parseOptions (String[] args) {
		CommandLine cmd = null;
		Options options = new Options();
		
		// Mandatory
		Option optionBam = Option.builder("b")
				.longOpt("bam").argName("bam file")
				.hasArg()
				.required(true)
				.desc("bam file path")
				.build();
		
		Option optionOutput= Option.builder("o")
				.longOpt("output").argName("string")
				.hasArg()
				.required(true)
				.desc("output file path")
				.build();
		
		
		
		options.addOption(optionBam)
		.addOption(optionOutput);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
	    ArrayList<String> tmpArgs = new ArrayList<String>();
	    for(int i=0; i<args.length; i++) {
	    	if(args[i].equalsIgnoreCase("-b") || args[i].equalsIgnoreCase("--bam") ||
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
		    
		    if(cmd.hasOption("b")) {
		    	bamFile = new File(cmd.getOptionValue("b"));
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
