package zhanglab.inputconvertor;

import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.systemsbiology.jrap.stax.MSXMLParser;
import org.systemsbiology.jrap.stax.Scan;

public class mzmlTest {

	public static void main(String[] args) throws IOException {
		Pattern fileNamePattern = Pattern.compile("[^\\/]+[$.]+(raw|mzml)", Pattern.CASE_INSENSITIVE);
		String fileName = "/Users/seunghyukchoi/eclipse-workspace/inputconvertor/VY20210311_FAIMS_HT_CPTAC1_W632_Fxn01.mzML";
		try {
			MSXMLParser parser = new MSXMLParser(fileName);
			int maxScanNum = parser.getMaxScanNumber();
			String originFilePath = parser.rapFileHeader().getParentFiles().get(0).getURI();
			Matcher matcher = fileNamePattern.matcher(originFilePath);
			matcher.find();
			String group = matcher.group();
			System.out.println(maxScanNum);
			
			for(int i=1; i<maxScanNum; i++) {
				Scan scan = parser.rap(i);
				int scanNum = scan.getHeader().getNum();
				int msLevel = scan.getHeader().getMsLevel();
				double ce = scan.getHeader().getCollisionEnergy();
				
				if(msLevel == 2) {
					System.out.println(scanNum+"\t"+ce);
					System.out.println(group);
				}
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
