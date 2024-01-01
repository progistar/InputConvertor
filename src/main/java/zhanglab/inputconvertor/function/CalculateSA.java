package zhanglab.inputconvertor.function;

import java.text.DecimalFormat;

import zhanglab.inputconvertor.data.Spectrum;

public class CalculateSA {
	
	public static double calculate (Spectrum s1, Spectrum s2, double tol, int precursorCharge, int chargeMax) {
		
		double[][] s1Peaks = s1.getAnnotatedPeaks(tol, precursorCharge, chargeMax);
		double[][] s2Peaks = s2.getAnnotatedPeaks(tol, precursorCharge, chargeMax);
		
		// normalization
		double s1PowerSum = getPowerSum(s1Peaks);
		double s2PowerSum = getPowerSum(s2Peaks);
		normalize(s1Peaks, Math.sqrt(s1PowerSum));
		normalize(s2Peaks, Math.sqrt(s2PowerSum));
		
		int s1Idx = 0;
		int s2Idx = 0;
		double ip = 0;
		boolean doMatchFurther = true;
		
		if(s1Peaks.length == 0 || s2Peaks.length == 0) {
			return 0;
		}
		
		while(doMatchFurther) {
			double delta = s1Peaks[s1Idx][0] - s2Peaks[s2Idx][0];
			if(Math.abs(delta) <= tol) {
				ip += s1Peaks[s1Idx][1] * s2Peaks[s2Idx][1];
				s1Idx++;
				s2Idx++;
			}
			else if(delta < 0) s1Idx++;
			else if(delta > 0) s2Idx++;
			
			if(s1Idx >= s1Peaks.length || s2Idx >= s2Peaks.length) doMatchFurther = false;
		}
		
		DecimalFormat decimalFormat = new DecimalFormat("#.#####");
		ip = Double.parseDouble(decimalFormat.format(ip));
		double sca = 1 - 2 * (Math.acos(ip) / Math.PI);
		
		return sca;
		
	}
	
	private static double getPowerSum (double[][] peaks) {
		double powerSum = 0;
		
		for(int i=0; i<peaks.length; i++) {
			powerSum += peaks[i][1] * peaks[i][1];
		}
		
		return powerSum;
	}
	
	private static void normalize (double[][] peaks, double normFactor) {
		for(int i=0; i<peaks.length; i++) {
			peaks[i][1] /= normFactor;
		}
	}

}
