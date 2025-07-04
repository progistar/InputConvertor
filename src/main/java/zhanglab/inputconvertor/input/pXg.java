package zhanglab.inputconvertor.input;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;

import zhanglab.inputconvertor.env.InputConvertorConstants;


public class pXg {
	public pXg () {}

	/**
	 * @deprecated
	 * 
	 * @param cmd
	 * @throws IOException
	 * @throws ParseException
	 */
	public void toNetMHCpanInputFormat(CommandLine cmd) throws IOException, ParseException {
		String inputFile = cmd.getOptionValue("i");
        String inputPattern = cmd.getOptionValue("p");

        int infPeptideIdx		= -1;
        
        File[] files = new File(inputFile).listFiles();
        if(files == null) {
        	System.out.println("-i option should be a path of folder including .fdr files.");
        	System.exit(1);
        }
        
        for(File file : files) {
        	if(file.getName().startsWith(".")) continue;
        	if(file.getName().endsWith(".fdr") && file.getName().contains(inputPattern)) {
        		
        		String outputName = inputFile+"/"+file.getName().replace(".fdr", ".netMHCpan.input");
        		BufferedWriter BW = new BufferedWriter(new FileWriter(outputName));
        		BufferedReader BR = new BufferedReader(new FileReader(file));
        		
        		String line = null;
        		String[] header = BR.readLine().split("\t");
        		
        		infPeptideIdx = InputConvertorConstants.getFieldIndex(header, InputConvertorConstants.IC_INFERRED_PEPTIDE_FIELD_NAME);
        		
        		Hashtable<String, String> isDuplicated = new Hashtable<String, String>();
        		while((line = BR.readLine()) != null) {
        			String[] fields = line.split("\t");
        			String infPeptide = fields[infPeptideIdx];
        			
        			if(infPeptide.length() >= 8 && infPeptide.length() <= 15) {

            			if(isDuplicated.get(infPeptide) == null) {
            				BW.append(infPeptide);
            				BW.newLine();
            				isDuplicated.put(infPeptide, "");
            			}
        			}
        			
        		}
        		
        		BR.close();
        		BW.close();
        		
        	}
        }
        
	}
	
}