package zhanglab.inputconvertor.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import org.systemsbiology.jrap.stax.MSXMLParser;
import org.systemsbiology.jrap.stax.Scan;

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
			SimpleMGFSelector(file);
		} else if(
				this.fileName.toLowerCase().endsWith(".mzml") || 
				this.fileName.toLowerCase().endsWith(".mzxml")) {
			SimpleMzMLSelector(file);
		}
	}
	
	public void SimpleMGFSelector (File file) throws IOException {
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
	
	public void SimpleMzMLSelector (File file) throws IOException {
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
