package zhanglab.inputconvertor.data;

public class HarmonyProtein implements Comparable<HarmonyProtein> {

	public String proteinId;
	public String geneName;
	public int proteinEvidence;
	
	public HarmonyProtein (String proteinId, String geneName, int proteinEvidence) {
		this.proteinId = proteinId;
		this.geneName = geneName;
		this.proteinEvidence = proteinEvidence;
	}

	@Override
	public int compareTo(HarmonyProtein o) {
		if(this.proteinEvidence > o.proteinEvidence) {
			return 1;
		} else if(this.proteinEvidence < o.proteinEvidence) {
			return -1;
		}
		return 0;
	}
}
