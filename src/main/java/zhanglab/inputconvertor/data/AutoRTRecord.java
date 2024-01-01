package zhanglab.inputconvertor.data;

public class AutoRTRecord implements Comparable<AutoRTRecord>{

	public double score;
	public String modifiedPeptide;
	public String rt;
	public String predRT;
	public int idx;
	public String fullRecord;
	
	@Override
	public int compareTo(AutoRTRecord o) {
		if(this.idx < o.idx) {
			return -1;
		}else if(this.idx > o.idx) {
			return 1;
		}
		return 0;
	}
	
}
