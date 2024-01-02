package zhanglab.inputconvertor.data;

import java.util.ArrayList;

public class NetMHCpanData {

	public String peptide;
	public ArrayList<HLA> hlas = new ArrayList<HLA>(); 
	
	public String getBestHLAType () {
		int bestHLAIndex = 0;
		
		for(int i=1; i<hlas.size(); i++) {
			if(hlas.get(i).elRank < hlas.get(bestHLAIndex).elRank) {
				bestHLAIndex = i;
			}
		}
		
		return hlas.get(bestHLAIndex).type;
	}
	
	public double getBestScore () {
		int bestHLAIndex = 0;
		
		for(int i=1; i<hlas.size(); i++) {
			if(hlas.get(i).elRank < hlas.get(bestHLAIndex).elRank) {
				bestHLAIndex = i;
			}
		}
		
		return hlas.get(bestHLAIndex).elRank;
	}
}
