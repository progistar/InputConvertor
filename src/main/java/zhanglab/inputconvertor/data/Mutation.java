package zhanglab.inputconvertor.data;

public class Mutation {

	public boolean isSomatic= false;
	public String refSeq	= null;
	public String altSeq	= null;
	public String chr		= null;
	public int pos			= -1;
	public byte type		= 0;
	
	
	public String toString() {
		return chr+":"+pos+"\t"+refSeq+">"+altSeq;
	}
}
