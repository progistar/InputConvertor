package zhanglab.inputconvertor.function;

import java.util.ArrayList;

import zhanglab.inputconvertor.env.InputConvertorConstants;

public abstract class Mode {
	
	public static ArrayList<Character> getStrandedness (int flags, String strandedness) {
		ArrayList<Character> strands = new ArrayList<Character>();
		boolean isFirstSegment = (0x40 & flags) == 0x40 ? true : false;
		boolean isForward = (0x10 & flags) == 0x10 ? false : true;
		
		// non-stranded
		if(strandedness.equalsIgnoreCase(InputConvertorConstants.NON_STRANDED)) {
			strands.add('+');
			strands.add('-');
		}  
		// Single-end
		else if(strandedness.equalsIgnoreCase(InputConvertorConstants.F_STRANDED)) {
			if(isForward) {
				strands.add('+');
			} else {
				strands.add('-');
			}
		} else if(strandedness.equalsIgnoreCase(InputConvertorConstants.R_STRANDED)) {
			if(isForward) {
				strands.add('-');
			} else {
				strands.add('+');
			}
		} 
		// Paired-end
		else {
			// R1
			if(isFirstSegment) {
				if(strandedness.equalsIgnoreCase(InputConvertorConstants.FR_STRANDED)) {
					if(isForward) {
						strands.add('+');
					} else {
						strands.add('-');
					}
				} else if(strandedness.equalsIgnoreCase(InputConvertorConstants.RF_STRANDED)) {
					if(isForward) {
						strands.add('-');
					} else {
						strands.add('+');
					}
				}
			} 
			// R2
			else {
				if(strandedness.equalsIgnoreCase(InputConvertorConstants.FR_STRANDED)) {
					if(isForward) {
						strands.add('-');
					} else {
						strands.add('+');
					}
				} else if(strandedness.equalsIgnoreCase(InputConvertorConstants.RF_STRANDED)) {
					if(isForward) {
						strands.add('+');
					} else {
						strands.add('-');
					}
				} 
			}
		}
		
		
		
		return strands;
	}
	
}
