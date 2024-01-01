package zhanglab.inputconvertor.module;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import org.apache.commons.cli.CommandLine;

import zhanglab.inputconvertor.data.SimpleMGFSelector;

public class SpectrumCount {

	public SpectrumCount (CommandLine cmd) throws IOException {
        String mgfFileBase = cmd.getOptionValue("f");
        
        File[] files = new File(mgfFileBase).listFiles();
        
        Hashtable<String, ArrayList<SimpleMGFSelector>> batch = new Hashtable<String, ArrayList<SimpleMGFSelector>>();
        for(File file : files) {
        	if(file.getName().startsWith(".")) continue;
        	if(file.getName().endsWith(".mgf")) {
        		
        		String[] fields = file.getName().split("\\_");
        		String key = fields[2]+"_"+fields[3];
        		
        		SimpleMGFSelector mgf = new SimpleMGFSelector(file);
        		
        		ArrayList<SimpleMGFSelector> single = batch.get(key);
        		if(single == null) {
        			single = new ArrayList<SimpleMGFSelector>();
        			batch.put(key, single);
        		}
        		single.add(mgf);
        	}
        }
        
        System.out.println("ic_identifier\tNo. files\tNo. spectra");
        batch.forEach((b, mgfList)->{
        	int cnt = 0;
        	for(SimpleMGFSelector mgf : mgfList) {
        		cnt += mgf.titles.size();
        	}
        	
        	System.out.println(b+"\t"+mgfList.size()+"\t"+cnt);
        	
        });
	}
}
