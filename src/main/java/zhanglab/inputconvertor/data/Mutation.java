package zhanglab.inputconvertor.data;

import zhanglab.inputconvertor.env.InputConvertorConstants;

public class Mutation {

	public boolean isSomatic= false;
	public String refSeq	= null;
	public String altSeq	= null;
	public String chr		= null;
	public int pos			= -1;
	public byte type		= 0;
	
	
	public void setMutationStatus () {
		if(refSeq.length() == altSeq.length()) {
			this.type = InputConvertorConstants.SNP;
		} else if(refSeq.length() > altSeq.length()) {
			this.type = InputConvertorConstants.DEL;
		} else if(refSeq.length() < altSeq.length()) {
			this.type = InputConvertorConstants.INS;
		}
	}
}
