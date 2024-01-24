package zhanglab.inputconvertor.data;

public class Mutation {

	public boolean isSomatic= false;
	public String refSeq	= null;
	public String altSeq	= null;
	public String chr		= null;
	public int pos			= -1;
	public byte type		= 0;
	public String key		= null;
	
	public String toString() {
		if(isSomatic) {
			return chr+":"+pos+"\t"+refSeq+">>"+altSeq;
		} else {
			return chr+":"+pos+"\t"+refSeq+">"+altSeq;
		}
	}
}
