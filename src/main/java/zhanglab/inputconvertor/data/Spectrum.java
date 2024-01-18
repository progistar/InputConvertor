package zhanglab.inputconvertor.data;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;

import zhanglab.inputconvertor.env.ProteomeConstants;

public class Spectrum {
	
	// 0: unsorted
	// 1: mz sorted, increasing order
	// -1: mz sorted, decreasing order
	// 2: int sorted, increasing order
	// -2: int sorted, decreasing order
	public int sortStatus = 0;

	public static final double INVALID_PEAK = Double.MAX_VALUE;
	
	public String title = null;
	public int scanNum = -1;
	public int charge = 0;
	public int msLevel = 0;
	public int index = 0; // this is one-based.
	public double retentionTime = 0;
	public double precursorMz = 0;
	public double precursorInt = 0;
	public double[] basePeak = null;
	public ArrayList<double[]> peaks = null;

	// ad-hoc
	public Peptide peptide = null;
	
	public Spectrum deepCopy() {
		ArrayList<double[]> newPeaks = new ArrayList<double[]>();
		Spectrum s = new Spectrum(scanNum, charge, msLevel, precursorMz, newPeaks, retentionTime, index);
		s.peptide = this.peptide;
		for(double[] peak : this.peaks) {
			double[] newPeak = {peak[0], peak[1]};
			newPeaks.add(newPeak);
		}
		
		return s;
	}
	
	/**
	 * If you want to skip putting peak list, then just give null.<br>
	 * 
	 * @param scanNum
	 * @param charge
	 * @param msLevel
	 * @param precursorMz
	 * @param peakList
	 * @param retentionTime
	 */
	public Spectrum (int scanNum, int charge, int msLevel, double precursorMz, ArrayList<double[]> peakList, double retentionTime, int index) {
		this.scanNum = scanNum;
		this.charge = charge;
		this.msLevel = msLevel;
		this.precursorMz = precursorMz;
		this.retentionTime = retentionTime;
		this.index = index;
		this.peaks = peakList;
		if(this.peaks == null) this.peaks = new ArrayList<double[]>();
	}
	
	/**
	 * Constructor of mzXML.
	 * 
	 * 
	 * @param scanNum
	 * @param charge
	 * @param msLevel
	 * @param precursorMz
	 * @param peakList
	 * @param retentionTime
	 * @param index
	 */
	public Spectrum (int scanNum, int charge, int msLevel, double precursorMz, double[][] peakList, double retentionTime, int index) {
		this.scanNum = scanNum;
		this.charge = charge;
		this.msLevel = msLevel;
		this.precursorMz = precursorMz;
		this.retentionTime = retentionTime;
		this.index = index;
		this.loadPeaksFromMzXML(peakList);
	}
	
	public double getMass () {
		double mass = this.charge * (this.precursorMz - ProteomeConstants.Proton);
		
		return mass;
	}

	public double getTotalIonChromatogram () {
		double TIC = .0;
		
		for(double[] peak : this.peaks) {
			TIC += peak[1];
		}
		
		return TIC;
	}
	
	
	/**
	 * loading peak information from mzXML.<br>
	 * 
	 * 
	 * 
	 * @param peaks
	 */
	public void loadPeaksFromMzXML (double[][] peaks) {
		if(this.peaks == null) this.peaks = new ArrayList<double[]>();
		int length = peaks[0].length;
		for(int i=0; i<length; i++) {
			double[] peak = new double[2];
			peak[0] = peaks[0][i];
			peak[1] = peaks[1][i];
			this.peaks.add(peak);
		}
	}
	/**
	 * mzOrInt: mz order = 0, int order = 1.<br>
	 * 
	 * @param mzOrInt
	 * @param isIncreasingOrder
	 */
	public void sortPeaks (final int mzOrInt, final boolean isIncreasingOrder) {
		int sortStatus = mzOrInt+1;
		if(!isIncreasingOrder) {
			sortStatus *= -1;
		}
		
		if(sortStatus == this.sortStatus) {
			return;
		} else {
			this.sortStatus = sortStatus;
		}
		
		Collections.sort(this.peaks, new Comparator<double[]>() {
			@Override
			public int compare(double[] peak1, double[] peak2) {
				if(peak1[mzOrInt] > peak2[mzOrInt]) return isIncreasingOrder ? 1 : -1;
				else if(peak1[mzOrInt] < peak2[mzOrInt]) return isIncreasingOrder ? -1 : 1;
				return 0;
			}
		});
	}
	
	public double[] getBasePeak () {
		int sizeOfPeaks = this.peaks.size();
		this.basePeak = this.peaks.get(0);
		for(int i=1; i<sizeOfPeaks; i++) {
			if(this.peaks.get(i)[1] > this.basePeak[1]) {
				this.basePeak = this.peaks.get(i);
			}
		}
		return this.basePeak;
	}
	
	/**
	 * get subpeaks using binary search.
	 * 
	 * @param fromMz
	 * @param toMz
	 * @return
	 */
	public ArrayList<double[]> getSubPeaks (double fromMz, double toMz) {
		ArrayList<double[]> subPeaks = new ArrayList<double[]>();
		this.sortPeaks(0, true);
		
		int sizeOfPeaks = this.peaks.size();
		int leftBound = 0;
		int rightBound = sizeOfPeaks-1;
		int startIndex = 0;
		while(leftBound <= rightBound) {
			int mid = (leftBound+rightBound) / 2;
			startIndex = mid;
			if(this.peaks.get(mid)[0] > fromMz) {
				rightBound = mid-1;
			} else if(this.peaks.get(mid)[0] < fromMz) {
				leftBound = mid+1;
			} else break;
		}
		
		startIndex -= 1;
		if(startIndex < 0) startIndex = 0;
		
		while(startIndex < sizeOfPeaks) {
			double mz = this.peaks.get(startIndex)[0];
			if(mz >= fromMz && mz <= toMz) {
				subPeaks.add(this.peaks.get(startIndex));
			} else if(mz > toMz){
				break;
			}
			startIndex++;
		}
		
		return subPeaks;
	}
	
	private double getSUMInt (ArrayList<Double> peaks) {
		double sum = 0;
		for(int i=0; i<peaks.size(); i++) {
			sum += peaks.get(i);
		}
		
		return sum;
	}
	public double[] getAnnotatedPeaks(double tolerance, int precursorCharge, int chargeLimit) {
		ArrayList<Double> annotatedPeaks = getPeaks(tolerance, precursorCharge, chargeLimit, true);
		ArrayList<Double> annotatedPeaksWithoutMods = getPeaks(tolerance, precursorCharge, chargeLimit, false);
		
		// Of course, it is right to consider modifications.
		// however, current version of Prosit does not reflect modification mass on each peak (I don't know why. It might be an error).
		// to cover those errors, I selected peaks having more intensity among two version of annotated peaks.
		
		double sumOfWMods = getSUMInt(annotatedPeaks);
		double sumOfWoMods = getSUMInt(annotatedPeaksWithoutMods);
		
		if(sumOfWMods < sumOfWoMods) {
			annotatedPeaks = annotatedPeaksWithoutMods;
		}
		
		double[] newPeaks = new double[annotatedPeaks.size()];
		for(int i=0; i<annotatedPeaks.size(); i++) {
			newPeaks[i] = annotatedPeaks.get(i);
		}
		
		return newPeaks;
	}
	
	private ArrayList<Double> getPeaks (double tolerance, int precursorCharge, int chargeLimit, boolean takingMods) {
		Hashtable<Double, Boolean> banDuplications = new Hashtable<Double, Boolean>();
		ArrayList<Double> annotatedPeaks = new ArrayList<Double>();

		// restricted to 1/2 ions.
		int upTo = precursorCharge-1;
		if(upTo == 0) {
			upTo = 1;
		} else if(upTo > chargeLimit) {
			upTo = chargeLimit;
		}

		
		for(int c=1; c<=upTo; c++) {
			double[] peaks = peptide.getTheoreticalLadder(ProteomeConstants.Y_ION, c, takingMods);
			for(int i=0; i<peaks.length; i++) {
				ArrayList<double[]> subPeaks = this.getSubPeaks(peaks[i]-tolerance, peaks[i]+tolerance);
				// select closest one
				double best = tolerance*2;
				double[] bestPeak = null;
				for(double[] peak : subPeaks) {
					if(Math.abs(peak[0]-peaks[i]) < best) {
						bestPeak = peak;
						best = Math.abs(peak[0]-peaks[i]);
					}
				}
				
				// if the peak is identified, then assign intensity.
				if(bestPeak != null && banDuplications.get(bestPeak[0]) == null) {
					banDuplications.put(bestPeak[0], true);
					annotatedPeaks.add(bestPeak[1]);
				} 
				// if there is no matched peak, than assign 0.
				else {
					annotatedPeaks.add(0.0);
				}
			}
			
			peaks = peptide.getTheoreticalLadder(ProteomeConstants.B_ION, c, takingMods);
			for(int i=0; i<peaks.length; i++) {
				ArrayList<double[]> subPeaks = this.getSubPeaks(peaks[i]-tolerance, peaks[i]+tolerance);
				// select closest one
				double best = tolerance*2;
				double[] bestPeak = null;
				for(double[] peak : subPeaks) {
					if(Math.abs(peak[0]-peaks[i]) < best) {
						bestPeak = peak;
						best = Math.abs(peak[0]-peaks[i]);
					}
				}
				
				// if the peak is identified, then assign intensity.
				if(bestPeak != null && banDuplications.get(bestPeak[0]) == null) {
					banDuplications.put(bestPeak[0], true);
					annotatedPeaks.add(bestPeak[1]);
				} 
				// if there is no matched peak, than assign 0.
				else {
					annotatedPeaks.add(0.0);
				}
			}
		}
		
		return annotatedPeaks;
	}
}