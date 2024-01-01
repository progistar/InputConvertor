package zhanglab.inputconvertor.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

public class SimpleMGFSelector {
	
	public String fileName = null;
	
	// (title, RT) pair
	public Hashtable<String, String> titleToRT = new Hashtable<String, String>();
	public ArrayList<String> titles = new ArrayList<String>();
	// (scanNum, title) pair
	public Hashtable<Integer, String> scanToTitle = new Hashtable<Integer, String>();
	
	public SimpleMGFSelector (File file) throws IOException {
		this.fileName = file.getName();
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		System.out.println("read: "+this.fileName);
		String line = null;
		
		String title = null;
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
				scanToTitle.put(Integer.parseInt(title.split("\\.")[len-2]), title);
			}
		}
		
		BR.close();
		
		System.out.println("A total of spectra: "+titleToRT.size());
	}
}
