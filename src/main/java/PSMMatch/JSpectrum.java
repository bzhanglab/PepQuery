package main.java.PSMMatch;

import java.util.ArrayList;
import java.util.Collections;


public class JSpectrum {

	public void setSpectrumTitle(String spectrumTitle) {
		this.spectrumTitle = spectrumTitle;
	}

	private String spectrumTitle;
	private String scanNumber;
	private double rt;
	private int charge = 0;



	private double parentMassToCharge;

	public double getParentMassToCharge() {
		return parentMassToCharge;
	}

	public void setParentMassToCharge(double parentMassToCharge) {
		this.parentMassToCharge = parentMassToCharge;
	}

	public double getParentMass() {
		return parentMass;
	}

	public void setParentMass(double parentMass) {
		this.parentMass = parentMass;
	}

	private double parentMass;
	private double intensity;
	public ArrayList<JPeak> peaks = new ArrayList<>();
	private ArrayList<JPeak> rawPeaks= new ArrayList<>();


	public JSpectrum(){
		
	}
	
	public String getSpectrumTitle() {
		return spectrumTitle;
	}

	public double getRt() {
		return rt;
	}
	public void setRt(double rt) {
		this.rt = rt;
	}
	public int getCharge() {
		return charge;
	}
	public void setCharge(int charge) {
		this.charge = charge;
	}

	public double getIntensity() {
		return intensity;
	}
	public void setIntensity(double intensity) {
		this.intensity = intensity;
	}
	public ArrayList<JPeak> getPeaks() {
		if(this.peaks.size()==0){
			this.resetPeaks();
		}
		return this.peaks;
	}

	public void setPeaks(ArrayList<JPeak> peaks) {
		this.peaks = peaks;
	}
	

	public void addRawPeak(JPeak peak){
		this.rawPeaks.add(peak);
	}


	public void resetPeaks(){
		this.peaks = new ArrayList<>();
		for(JPeak jPeak :this.rawPeaks){
			JPeak tmpJPeak = jPeak.clone();
			this.peaks.add(tmpJPeak);
		}
		
	}
	

	public void sortPeaksByMZ(){
		PeakMZComparator pComparator = new PeakMZComparator();
		Collections.sort(this.peaks, pComparator);
	}

	public void sortPeaksByIntensity(){
		PeakIntensityComparator pComparator = new PeakIntensityComparator();
		Collections.sort(this.peaks,pComparator);
	}

	public double getMaxIntensity(){
		double maxIntensity= 0;
		for(JPeak jPeak:this.peaks){
			if(maxIntensity<jPeak.getIntensity()){
				maxIntensity = jPeak.getIntensity();
			}
		}
		return maxIntensity;
	}


	public double getTotalIonCurrent(){
		double totalIonCurrent = 0.0;
		for(JPeak jPeak: this.peaks){
			totalIonCurrent+=jPeak.getIntensity();
		}
		return totalIonCurrent;
	}


	
}
