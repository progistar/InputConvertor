package zhanglab.inputconvertor.function;

import java.text.DecimalFormat;

import zhanglab.inputconvertor.data.Spectrum;

public class CalculateSA {
	
	public static double calculate (Spectrum s1, Spectrum s2, double tol, int precursorCharge, int chargeMax) {
		
		double[] s1Peaks = s1.getAnnotatedPeaks(tol, precursorCharge, chargeMax);
		double[] s2Peaks = s2.getAnnotatedPeaks(tol, precursorCharge, chargeMax);
		
		// normalization
		double s1PowerSum = getPowerSum(s1Peaks);
		double s2PowerSum = getPowerSum(s2Peaks);
		normalize(s1Peaks, Math.sqrt(s1PowerSum));
		normalize(s2Peaks, Math.sqrt(s2PowerSum));
		
		double ip = 0;
		if(s1PowerSum == 0 || s2PowerSum == 0) {
			return 0;
		}
		
		for(int i=0; i<s1Peaks.length; i++) {
			ip += s1Peaks[i] * s2Peaks[i];
		}
		
		DecimalFormat decimalFormat = new DecimalFormat("#.#####");
		ip = Double.parseDouble(decimalFormat.format(ip));
		double sca = 1 - 2 * (Math.acos(ip) / Math.PI);
		
		return sca;
		
	}
	
	private static double getPowerSum (double[] peaks) {
		double powerSum = 0;
		
		for(int i=0; i<peaks.length; i++) {
			powerSum += peaks[i] * peaks[i];
		}
		
		return powerSum;
	}
	
	private static void normalize (double[] peaks, double normFactor) {
		for(int i=0; i<peaks.length; i++) {
			peaks[i] /= normFactor;
		}
	}

}
