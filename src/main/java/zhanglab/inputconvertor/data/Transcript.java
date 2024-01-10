package zhanglab.inputconvertor.data;

import java.util.ArrayList;

public class Transcript {

	public String tID;
	public String strand;
	public String chr;
	public String start;
	public String end;
	public String attrs;
	
	public ArrayList<Exon> exons = new ArrayList<Exon>();
}
