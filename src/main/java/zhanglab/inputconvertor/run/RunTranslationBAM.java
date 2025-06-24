package zhanglab.inputconvertor.run;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
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
	public static File logFile = null;
	public static int maxMer = 25;
	public static double rphm = 10;
	public static int atMost = 100000;
	
	public static void main(String[] args) throws IOException {
		parseOptions(args);
		
		System.out.println("RPHM threshold: "+rphm);
		System.out.println("Minimum peptide length: "+InputConvertorConstants.MIN_PEPT_LEN);
		System.out.println("Maximum peptide length: "+maxMer);
		System.out.println("Maximum number of entries: "+atMost);
		System.out.println("## Unmapped translation ##");
		System.out.println("Load "+bamFile.getName());
		writeFastaEntry();
	}
	
	public static void writeFastaEntry () throws IOException {
		ArrayList<FastaEntry> fastaEntries = new ArrayList<FastaEntry>();
		BufferedWriter log = new BufferedWriter(new FileWriter(logFile));
		String strandedness = "NON";
		
		int R1F = 0;
		int R1R = 0;
		int R2F = 0;
		int R2R = 0;
		int totalReads = 0;
		
		try (SamReader samReader = SamReaderFactory.makeDefault().open(bamFile)) {
			SAMRecordIterator iterator = samReader.iterator();
			while(iterator.hasNext()) {
				SAMRecord samRecord = iterator.next();
				
				boolean isPass = false;
				if(samRecord.isSecondaryAlignment()) {
            		isPass = true;
            	} else {
            		totalReads ++;
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
			
			
			System.out.println("## Total primary reads: "+totalReads);
			
			log.append("## Total primary reads: "+totalReads);
			log.newLine();
			log.append("## 1F:1R:2F:2R="+R1F+":"+R1R+":"+R2F+":"+R2R);
			log.newLine();
			log.append("## "+strandedness);
			log.newLine();
			
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		int totalUmap = 0;
		try (SamReader samReader = SamReaderFactory.makeDefault().open(bamFile)) {
			SAMRecordIterator iterator = samReader.queryUnmapped();
			while (iterator.hasNext()) {
				iterator.next();
				totalUmap++;
			}
			System.out.println("Total unmapped reads: "+totalUmap +" (" + 100.0*(totalUmap+0.0)/(totalReads)+"%)");
			iterator.close();
			
			log.append("## Total unmapped reads: "+totalUmap +" (" + 100.0*(totalUmap+0.0)/(totalReads)+"%)");
			log.newLine();
		}
		
		try (SamReader samReader = SamReaderFactory.makeDefault().open(bamFile)) {
			SAMRecordIterator iterator = samReader.queryUnmapped();
			Hashtable<String, Integer> readMap = new Hashtable<String, Integer>();
			
			double percent = 0;
			double cnt = 0;
			
			while (iterator.hasNext()) {
	            SAMRecord samRecord = iterator.next();
	            int flags = samRecord.getFlags();
	            ArrayList<Character> strands = Mode.getStrandedness(flags, strandedness);
	            
	            cnt++;
	            if(Math.abs(cnt/totalUmap - percent) < 0.00001) {
	            	System.out.println((int)(100*percent)+"%");
	            	percent+=0.01;
	            }
	            
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
			
			int[] readCounts = new int[1001];
			double[] readToRPHM = new double[1001];
			readCounts[0] = totalReads;
			readMap.forEach((pept, reads)->{
				int idx = Math.min(reads, readCounts.length-1);
				readCounts[idx]++;
				double thisRPHM = Math.pow(10, 8)*(idx+0.0)/(readCounts[0]);
				readToRPHM[idx] = thisRPHM;
				if(thisRPHM > rphm) {
					FastaEntry entry = new FastaEntry();
					entry.tool = InputConvertorConstants.UNMAPPED_HEADER_ID;
					entry.idx = fastaEntries.size()+1;
					entry.sequence = pept;
					entry.strand = ".";
					entry.frame = ".";
					entry.transcriptId = "Unmapped";
					entry.geneId = "" + (Math.pow(10, 8)*(reads+0.0)/(readCounts[0]));
					entry.geneName = ".";
					entry.description = ".";
					fastaEntries.add(entry);
				}
			});
			
			log.append("## Total enumerated peptides: "+readMap.size());
			log.newLine();
			log.append("## Peptides with RPHM > "+rphm+": "+fastaEntries.size()+" ("+100.0*(fastaEntries.size()+0.0)/(readMap.size())+"%)");
			log.newLine();
			
			log.append("Read\tRPHM\tCount");
			log.newLine();
			for(int i=1; i<readCounts.length; i++) {
				log.append(i+"\t"+readToRPHM[i]+"\t"+readCounts[i]);
				log.newLine();
			}
			
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		// sort by geneId (=RPHM)
		Collections.sort(fastaEntries, new Comparator<FastaEntry>() {
			@Override
			public int compare(FastaEntry o1, FastaEntry o2) {
				double rphm1 = Double.parseDouble(o1.geneId);
				double rphm2 = Double.parseDouble(o2.geneId);
						
				if(rphm1 > rphm2) {
					return -1;
				} else if(rphm1 < rphm2) {
					return 1;
				}
						
				return 0;
			}
		});
		
		
		BufferedWriter BW = new BufferedWriter(new FileWriter(outputFile));
	        
        for(FastaEntry entry : fastaEntries) {
        	BW.append(">"+entry.toMeta(null));
        	BW.newLine();
        	BW.append(entry.sequence);
        	BW.newLine();
        }
        
        BW.close();
        log.close();
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
		
		Option optionRPHM = Option.builder("r")
				.longOpt("rphm").argName("float")
				.hasArg()
				.required(false)
				.desc("RPHM threshold to print (Default is 10).")
				.build();
		
		Option optionAtMost = Option.builder("n")
				.longOpt("num").argName("natural number")
				.hasArg()
				.required(false)
				.desc("maximum number of entires to print (Default is 100,000).")
				.build();
		
		Option optionOutput= Option.builder("o")
				.longOpt("output").argName("string")
				.hasArg()
				.required(true)
				.desc("output file path")
				.build();
		
		
		
		options.addOption(optionBam)
		.addOption(optionRPHM)
		.addOption(optionAtMost)
		.addOption(optionOutput);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
	    ArrayList<String> tmpArgs = new ArrayList<String>();
	    for(int i=0; i<args.length; i++) {
	    	if(args[i].equalsIgnoreCase("-b") || args[i].equalsIgnoreCase("--bam") ||
			args[i].equalsIgnoreCase("-o") || args[i].equalsIgnoreCase("--output") ||
			args[i].equalsIgnoreCase("-n") || args[i].equalsIgnoreCase("--num") ||
			args[i].equalsIgnoreCase("-r") || args[i].equalsIgnoreCase("--rphm")) {
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
		    	logFile = new File(cmd.getOptionValue("o")+".log");
		    }
		    if(cmd.hasOption("r")) {
		    	rphm = Double.parseDouble(cmd.getOptionValue("r"));
		    }
		    if(cmd.hasOption("n")) {
		    	atMost = Integer.parseInt(cmd.getOptionValue("n"));
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
