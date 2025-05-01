package zhanglab.inputconvertor.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.systemsbiology.jrap.stax.MSXMLParser;
import org.systemsbiology.jrap.stax.Scan;

import zhanglab.inputconvertor.env.InputConvertorConstants;

public class Spectra {

	public ArrayList<Spectrum> spectra = new ArrayList<Spectrum>();
	// indexed by title
	public Hashtable<String, Spectrum> scanIndexer = new Hashtable<String, Spectrum>();
	//
	public boolean isComplete = false;
	public String predictType = InputConvertorConstants.EXPERIMENTAL_SPECTRUM;
	
	/**
	 * If the predicted spectra,
	 * 
	 * 
	 * 
	 * @param fileName
	 * @param isPredicted
	 */
	public Spectra (String fileName, String predictType) {
		this.predictType = predictType;
		
		
		if(fileName.toLowerCase().endsWith(".mgf")) {
			MGF(fileName);
		} else if(
				fileName.toLowerCase().endsWith(".mzml") || 
				fileName.toLowerCase().endsWith(".mzxml")) {
			MzML(fileName);
		}
	}
	
	private void MGF (String fileName) {
		try {
			BufferedReader BR = new BufferedReader(new FileReader(fileName));
			String line = null;
			int index = -1;
			double pepMass = 0;
			double precursorInt = 0;
			int charge = 0;
			double rt = 0;
			String title = null;
			int scanNum = -1;
			ArrayList<double[]> peaks = null;
			
			Pattern peakPattern = Pattern.compile("^[0-9]");
			Pattern scanPattern = Pattern.compile("[\\.]+[0-9]+[\\.]+");
			
			while((line = BR.readLine()) != null) {
				if(line.length() == 0) continue;
				
				if(line.startsWith("BEGIN")) {
					peaks = new ArrayList<double[]>();
					index ++;
				}else if(line.startsWith("TITLE")) {
					title = line.split("\\=")[1].split("\\s")[0];

					Matcher matcher = scanPattern.matcher(title);
					if(matcher.find()) {
						scanNum = Integer.parseInt(matcher.group().replace(".", ""));
					}
				}else if(line.startsWith("RTIN")) {
					rt = Double.parseDouble(line.split("\\=")[1]);
				}else if(line.startsWith("PEPMASS")) {
					pepMass = Double.parseDouble(line.split("\\s")[0].split("\\=")[1]);
					try {
						precursorInt = Double.parseDouble(line.split("\\s")[1]);
					}catch(Exception e) {
						precursorInt = -1;
					}
				}else if(line.startsWith("CHARGE")){
					charge = (int) Double.parseDouble(line.split("\\=")[1].replace("+", ""));
					if(charge == 0) {
						System.out.println(Double.parseDouble(line.split("\\=")[1].replace("+", "")));
					}
				}else if(peakPattern.matcher(line).find()) {
					double[] peak = new double[2];
					String[] peakStr = line.split("\\s");
					peak[0] = Double.parseDouble(peakStr[0]);
					peak[1] = Double.parseDouble(peakStr[1]);
					peaks.add(peak);
				}else if(line.startsWith("END")) {
					Spectrum spectrum = new Spectrum(scanNum, charge, 2, pepMass, peaks, rt, index);
					spectrum.precursorInt = precursorInt;
					spectrum.title = title;
					spectra.add(spectrum);
					
					// if scanNum is -1, then it is missing value.
					if(title != null) {
						scanIndexer.put(title, spectrum);
					}
				}
			}
			
			BR.close();
			isComplete = true;
		}catch(IOException ioe) {
			ioe.printStackTrace();
		}
	}
	
	private void MzML (String fileName) {
		try {
			MSXMLParser parser = new MSXMLParser(fileName);
			File file = new File(fileName);
			String titlePrefix = file.getName();
			titlePrefix = titlePrefix.substring(0, titlePrefix.lastIndexOf("."));
			int maxScanNum = parser.getMaxScanNumber();
			
			for (int i =1; i<=maxScanNum; i++) {
				Scan scan = parser.rap(i);
				if( scan == null) continue;
				double[][] peaks = scan.getMassIntensityList();
				int scanNum = scan.getHeader().getNum();
				int msLevel = scan.getHeader().getMsLevel();
				double retentionTime = scan.getHeader().getDoubleRetentionTime();
				
				if(msLevel == 1) {
					
				} 
				// MS/MS
				else if(msLevel == 2) {
					int precursorCharge = scan.getHeader().getPrecursorCharge();
					double precursorMz = scan.getHeader().getPrecursorMz();
					double precursorInt = scan.getHeader().getPrecursorIntensity();
					
					Spectrum spectrum = new Spectrum(scanNum, precursorCharge, msLevel, precursorMz, peaks, retentionTime, i);
					spectrum.title = titlePrefix+"."+scanNum+"."+scanNum+"."+precursorCharge;
					spectrum.precursorMz = precursorMz;
					spectrum.precursorInt = precursorInt;
					spectra.add(spectrum);
					
					// if scanNum is -1, then it is missing value.
					if(spectrum.title != null) {
						scanIndexer.put(spectrum.title, spectrum);
					}
				}
			}
			
			isComplete = true;
			
		} catch (Exception e) {
			
		}
	}
}