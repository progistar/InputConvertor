package zhanglab.inputconvertor.data;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import org.systemsbiology.jrap.stax.MSXMLParser;
import org.systemsbiology.jrap.stax.Scan;

import zhanglab.inputconvertor.env.InputConvertorConstants;
import zhanglab.inputconvertor.run.Run;

public class SimpleSpectraSelector {
	
	public String fileName = null;
	
	// (title, RT) pair
	public Hashtable<String, String> titleToRT = new Hashtable<String, String>();
	public ArrayList<String> titles = new ArrayList<String>();
	// (scanNum, title) pair
	public Hashtable<Integer, String> scanToTitle = new Hashtable<Integer, String>();
	// (index, scanNum) pair
	public Hashtable<Integer, Integer> indexToScan = new Hashtable<Integer, Integer>();
	
	public SimpleSpectraSelector (File file) throws IOException {
		this.fileName = file.getName();
		
		if(this.fileName.toLowerCase().endsWith(".mgf")) {
			if(Run.tool.equalsIgnoreCase(InputConvertorConstants.PEAKS12)) {
				simplePEAKS12MGFSelector(file);
			} else {
				simpleMGFSelector(file);
			}
			
		} else if(
				this.fileName.toLowerCase().endsWith(".mzml") || 
				this.fileName.toLowerCase().endsWith(".mzxml")) {
			simpleMzMLSelector(file);
		}
	}
	
	public void simpleMGFSelector (File file) throws IOException {
		this.fileName = file.getName();
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		System.out.println("read: "+this.fileName);
		String line = null;
		
		String title = null;
		int idx = 0;
		while((line = BR.readLine()) != null) {
			if(line.startsWith("TITLE=")) {
				title = line.split("\\=")[1].split("\\s")[0];
				if(titleToRT.get(title) != null) {
					System.out.println("warning: "+title+" is duplicated!");
				}
			} else if(line.startsWith("RTIN")) {
				String rt = line.split("\\=")[1];
				titleToRT.put(title, rt);
				titles.add(title);
				
				int len = title.split("\\.").length;
				int scanNum = Integer.parseInt(title.split("\\.")[len-2]);
				scanToTitle.put(scanNum, title);
				indexToScan.put(idx++, scanNum);
			}
		}
		
		BR.close();
		
		System.out.println("A total of spectra: "+titleToRT.size());
	}
	
	/**
	 * 
	 * For mgf files exported by PEAKS12,
	 * It is difficult to match results to spectrum using title.
	 * Instead, it is better to use a combination of precursor id and pepmass as a key.
	 * 
	 * 
	 * 
	 * @param file
	 * @throws IOException
	 */
	public void simplePEAKS12MGFSelector (File file) throws IOException {
		File outFile = new File(file.getAbsolutePath().replace(".mgf", ".peaks12.mgf"));
		this.fileName = file.getName();
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		System.out.println("read: "+this.fileName);
		String line = null;
		
		String precursorid = null;
		int mass = 0;
		ArrayList<String> spectrum = new ArrayList<String>();
		String newTitle = null;
		BufferedWriter BW = new BufferedWriter(new FileWriter(outFile));
		while((line = BR.readLine()) != null) {
			if(line.length() == 0) continue;
			
			if(line.startsWith("BEGIN")) {
				if(spectrum.size() > 0) {
					// write file
					// if title? replace
					for(String record : spectrum) {
						if(record.startsWith("TITLE=")) {
							record = "TITLE="+newTitle;
						}
						BW.append(record);
						BW.newLine();
					}
				}
				spectrum.clear();
			}
			
			spectrum.add(line);
			
			if(line.startsWith("PRECURSORID=")) {
				precursorid = line.split("\\=")[1].split("\\s")[0];
			}
			else if(line.startsWith("PEPMASS=")) {
				mass = (int) (1000 * Double.parseDouble(line.split("\\=")[1].split("\\s")[0]));
			}
			else if(line.startsWith("RTIN")) {
				newTitle = precursorid+"+"+mass;
				String rt = line.split("\\=")[1];
				titleToRT.put(newTitle, rt);
				titles.add(newTitle);
			}
		}
		
		if(spectrum.size() > 0) {
			// write file
			// if title? replace
			for(String record : spectrum) {
				if(record.startsWith("TITLE=")) {
					record = "TITLE="+newTitle;
				}
				BW.append(record);
				BW.newLine();
			}
		}
		spectrum.clear();
		
		BR.close();
		BW.close();
		
		System.out.println("A total of spectra: "+titleToRT.size());
	}
	
	public void simpleMzMLSelector (File file) throws IOException {
		this.fileName = file.getName();
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		System.out.println("read: "+this.fileName);
		String line = null;
		
		try {
			MSXMLParser parser = new MSXMLParser(file.getAbsolutePath());
			String titlePrefix = file.getName();
			titlePrefix = titlePrefix.substring(0, titlePrefix.lastIndexOf("."));
			int maxScanNum = parser.getMaxScanNumber();
			
			for (int i =1; i<=maxScanNum; i++) {
				Scan scan = parser.rap(i);
				if( scan == null) continue;
				int scanNum = scan.getHeader().getNum();
				int msLevel = scan.getHeader().getMsLevel();
				double retentionTime = scan.getHeader().getDoubleRetentionTime();
				
				if(msLevel == 1) {
					
				} 
				// MS/MS
				else if(msLevel == 2) {
					int precursorCharge = scan.getHeader().getPrecursorCharge();
					
					String title = titlePrefix+"."+scanNum+"."+scanNum+"."+precursorCharge;;
					if(titleToRT.get(title) != null) {
						System.out.println("warning: "+title+" is duplicated!");
					}
					
					String rt = retentionTime+"";
					titleToRT.put(title, rt);
					titles.add(title);
					scanToTitle.put(scanNum, title);
				}
			}
		} catch (Exception e) {
			
		}
		
		BR.close();
		
		System.out.println("A total of spectra: "+titleToRT.size());
	}
}
