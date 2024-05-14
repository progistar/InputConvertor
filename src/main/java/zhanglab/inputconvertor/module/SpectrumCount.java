package zhanglab.inputconvertor.module;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import org.apache.commons.cli.CommandLine;

import zhanglab.inputconvertor.data.SimpleSpectraSelector;

public class SpectrumCount {

	public SpectrumCount (CommandLine cmd) throws IOException {
        String mgfFileBase = cmd.getOptionValue("f");
        
        File[] files = new File(mgfFileBase).listFiles();
        
        Hashtable<String, ArrayList<SimpleSpectraSelector>> batch = new Hashtable<String, ArrayList<SimpleSpectraSelector>>();
        for(File file : files) {
        	if(file.getName().startsWith(".")) continue;
        	if(file.getName().endsWith(".mgf")) {
        		
        		String[] fields = file.getName().split("\\_");
        		String key = fields[2]+"_"+fields[3];
        		
        		SimpleSpectraSelector mgf = new SimpleSpectraSelector(file);
        		
        		ArrayList<SimpleSpectraSelector> single = batch.get(key);
        		if(single == null) {
        			single = new ArrayList<SimpleSpectraSelector>();
        			batch.put(key, single);
        		}
        		single.add(mgf);
        	}
        }
        
        System.out.println("ic_identifier\tNo. files\tNo. spectra");
        batch.forEach((b, mgfList)->{
        	int cnt = 0;
        	for(SimpleSpectraSelector mgf : mgfList) {
        		cnt += mgf.titles.size();
        	}
        	
        	System.out.println(b+"\t"+mgfList.size()+"\t"+cnt);
        	
        });
	}
}
