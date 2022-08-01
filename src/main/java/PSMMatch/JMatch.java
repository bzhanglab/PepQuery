package main.java.PSMMatch;

import com.compomics.util.experiment.biology.ions.Ion;
import com.compomics.util.experiment.biology.ions.NeutralLoss;
import com.compomics.util.experiment.biology.ions.impl.PeptideFragmentIon;
import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.identification.matches.IonMatch;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.experiment.identification.spectrum_annotation.AnnotationParameters;
import com.compomics.util.experiment.identification.spectrum_annotation.SpecificAnnotationParameters;
import com.compomics.util.experiment.identification.spectrum_annotation.SpectrumAnnotator;
import com.compomics.util.experiment.identification.spectrum_annotation.spectrum_annotators.PeptideSpectrumAnnotator;
import com.compomics.util.experiment.identification.spectrum_assumptions.PeptideAssumption;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;


import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import com.compomics.util.parameters.identification.search.ModificationParameters;
import main.java.pg.CParameter;
import main.java.pg.JSequenceProvider;
import main.java.pg.PeptideSearchMT;


import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.*;


/**
 * Scoring algorithm: basic functions used by different algorithms.
 *
 */
public class JMatch {


	/**
	 * The maximum number of peaks in a spectrum
	 */
	public int iMaxPeaks=50;
	/**
	 * The lowest m/z in a spectrum, a peak with m/z < dLowestMass will be removed from scoring step
	 */
	public double dLowestMass=150.0;

	/**
	 * The intensity cutoff, default is 0.05
	 */
	public double intensityLimit = 0.05;

	/**
	 * Used by MVH
	 */
	public double ticCutoffPercentage=0.95;


	public Peptide objPeptide = null;


	public Spectrum spectrum;
	public ArrayList<IonMatch> matches = new ArrayList<>();
	public JSpectrum jSpectrum = new JSpectrum();
	public JPeptideSpectrumMatch jPeptideSpectrumMatch = new JPeptideSpectrumMatch();

	/**
	 * peptide fragmentation method: cid/hcd
	 */
	public String fragmentMethod = "cid";

	/**
	 * The max charge state for fragment ion matching.
	 */
	public int maxFragmentChargeState=3;



	private boolean lossWaterNH3=false;

	public static double isotope_tol = 0.005;


	public JMatch(){
		
	}
	

	public JMatch(Spectrum spectrum, JSpectrum jSpectrum, JPeptideSpectrumMatch psm){
		this.spectrum = spectrum;
		this.jSpectrum = jSpectrum;
		this.jPeptideSpectrumMatch = psm;
	}

	public void initialize(boolean lossWaterNH3, int maxFragCharge) throws IllegalArgumentException, FileNotFoundException, ClassNotFoundException, IOException, InterruptedException, SQLException{
		setLossWaterNH3(lossWaterNH3);
		setMaxFragmentChargeState(maxFragCharge);
		// Peptide spectrum matching
		getSpectrumAnnotation();
		
	}

	public int removeParentPeak(){
		//mz for ms2 tol : unit is Da.
		int precursor_charge = 0;
		if(this.jPeptideSpectrumMatch != null){
			precursor_charge = this.jPeptideSpectrumMatch.getCharge();
		}else{
			precursor_charge = jSpectrum.getCharge();
		}

		if(precursor_charge ==0 && spectrum != null){
			precursor_charge = spectrum.getPrecursor().possibleCharges[0];
		}
		if(precursor_charge == 0){
			precursor_charge = 1;
		}
		double tol_mz = 1.0*CParameter.itol/precursor_charge;
		int rpn = 0;
		Iterator<JPeak> pListIterator = this.jSpectrum.getPeaks().iterator();  
		while(pListIterator.hasNext()){
			JPeak jPeak = pListIterator.next();
			if(Math.abs( jPeak.getMz() - this.jPeptideSpectrumMatch.getExperimentalMassToCharge() ) <= tol_mz){
				rpn++;
				pListIterator.remove();
			}
		}
		return rpn;
	}
	

	public int removeLowMasses(){
		int rpn = 0;
		Iterator<JPeak> pListIterator = this.jSpectrum.getPeaks().iterator();  
		while(pListIterator.hasNext()){
			JPeak jPeak = pListIterator.next();
			if( jPeak.getMz() <= this.dLowestMass){
				rpn++;
				pListIterator.remove();
			}
		}
		
		return rpn;
	}
	
	/**
	 * Keep the top iMaxPeaks fragment ions.
	 * @return
	 */
	public boolean filterByPeakCount(){
		boolean isRemove = false;
		int peaksNumber = this.jSpectrum.getPeaks().size();
		if(peaksNumber>this.iMaxPeaks){
			this.jSpectrum.sortPeaksByIntensity();
			this.jSpectrum.getPeaks().subList(0, peaksNumber-this.iMaxPeaks).clear();
			isRemove=true;
		}
		return isRemove;
	}
	
	/**
	 * Keep the top iMaxPeaks fragment ions.
	 * @return
	 */
	public boolean filterByPeakCount(int maxPeakCount){
		boolean isRemove = false;
		int peaksNumber = this.jSpectrum.getPeaks().size();
		if(peaksNumber>maxPeakCount){
			this.jSpectrum.sortPeaksByIntensity();
			this.jSpectrum.getPeaks().subList(0, peaksNumber-maxPeakCount).clear();
			isRemove=true;
		}
		return isRemove;
	}
	
	/**
	 * Remove low intensity fragment ions
	 * @return
	 */
	public int removeLowIntensityPeak(){
		int rpn = 0;
		double maxIntensity = this.jSpectrum.getMaxIntensity();
		double limitIntensity = 1.0*this.intensityLimit*maxIntensity;

		Iterator<JPeak> pListIterator = this.jSpectrum.getPeaks().iterator();  
		while(pListIterator.hasNext()){
			JPeak jPeak = pListIterator.next();
			if( jPeak.getIntensity()<limitIntensity){
				rpn++;
				pListIterator.remove();
			}
		}
		
		return rpn;
	}
	
	/**
	 * Remove low intensity fragment ions
	 * @param interLimit intensity cutoff
	 * @return the number of removed peaks
	 */
	public int removeLowIntensityPeak( double interLimit){
		int rpn = 0;
		double maxIntensity = this.jSpectrum.getMaxIntensity();
		double limitIntensity = 1.0*interLimit*maxIntensity;

		Iterator<JPeak> pListIterator = this.jSpectrum.getPeaks().iterator();  
		while(pListIterator.hasNext()){
			JPeak jPeak = pListIterator.next();
			if( jPeak.getIntensity()<limitIntensity){
				rpn++;
				pListIterator.remove();
			}
		}
		
		return rpn;
	}
	
	/**
	 * Remove isotopes removes multiple entries within 0.95 Da of each other, retaining
	 * the highest value. this is necessary because of the behavior of some peak 
	 * finding routines in commercial software
	 * @return the number of removed peaks.
	 */
	public int removeIsotopes(){
		int rpn = 0;
		if(this.jSpectrum.getPeaks().size()>2){
			// sort peaks by mz: from small to large
			this.getJSpectrum().sortPeaksByMZ();
			JPeak jPeak1 = this.getJSpectrum().getPeaks().get(0);
			JPeak jPeak2 = this.getJSpectrum().getPeaks().get(1);
			ArrayList< JPeak> tmpJPeaks = new ArrayList<JPeak>();
			double mz1 = jPeak1.getMz();
			for(int i=1;i<this.getJSpectrum().getPeaks().size();i++){
				jPeak2 = this.getJSpectrum().getPeaks().get(i);
				if(jPeak2.getMz() - mz1 >= 0.95 || jPeak2.getMz() < 200){
					tmpJPeaks.add(jPeak1);
					jPeak1 = jPeak2;
					mz1 = jPeak1.getMz();
				}else if (jPeak2.getIntensity() > jPeak1.getIntensity()){
					jPeak1 = jPeak2;
					mz1 = jPeak1.getMz();
				}
			}
			rpn = this.jSpectrum.getPeaks().size() - tmpJPeaks.size();
			tmpJPeaks.add(jPeak1);
			this.jSpectrum.setPeaks(tmpJPeaks);
		}
		return rpn;
	}
	
	
	/**
	 * clean_isotopes removes peaks that are probably C13 isotopes
	 * @return the number of removed peaks.
	 */
	public int cleanIsotopes(){
		int rpn = 0;
		if(this.jSpectrum.getPeaks().size()>2){
			// sort peaks by mz: from small to large
			this.getJSpectrum().sortPeaksByMZ();
			JPeak jPeak1 = this.getJSpectrum().getPeaks().get(0);
			JPeak jPeak2 = this.getJSpectrum().getPeaks().get(1);
			ArrayList< JPeak> tmpJPeaks = new ArrayList<JPeak>();
			double mz1 = jPeak1.getMz();
			for(int i=1;i<this.getJSpectrum().getPeaks().size();i++){
				jPeak2 = this.getJSpectrum().getPeaks().get(i);
				if(jPeak2.getMz() - mz1 >= 1.5 || jPeak2.getMz() < 200){
					tmpJPeaks.add(jPeak1);
					jPeak1 = jPeak2;
					mz1 = jPeak1.getMz();
				}else if (jPeak2.getIntensity() > jPeak1.getIntensity()){
					jPeak1 = jPeak2;
					mz1 = jPeak1.getMz();
				}
			}
			rpn = this.jSpectrum.getPeaks().size() - tmpJPeaks.size();
			tmpJPeaks.add(jPeak1);
			this.jSpectrum.setPeaks(tmpJPeaks);
		}
		return rpn;
	}
	
	/**
	 * Normalize fragment ion intensities.
	 * @param maxPeakAfterNornalize The max peak intensity after normalization
	 */
	public void normalizePeaks(double maxPeakAfterNornalize){
		double maxIntensity = this.jSpectrum.getMaxIntensity();
		for(int i=0;i<this.jSpectrum.getPeaks().size();i++){
			double norIntensity = 1.0*maxPeakAfterNornalize*this.jSpectrum.getPeaks().get(i).getIntensity()/maxIntensity;
			this.jSpectrum.getPeaks().get(i).setIntensity(norIntensity);
		}
	}

	public Spectrum getMSnSpectrum() {
		return spectrum;
	}


	public void setMSnSpectrum(Spectrum spectrum) {
		this.spectrum = spectrum;
	}
	

	public ArrayList<IonMatch> getIonMatches() {
		return matches;
	}


	public JSpectrum getJSpectrum() {
		return jSpectrum;
	}
	

	public void setJSpectrum(JSpectrum jSpectrum) {
		this.jSpectrum = jSpectrum;
	}


	public String getFragmentMethod() {
		return fragmentMethod;
	}

	public void setFragmentMethod(String fragmentMethod) {
		this.fragmentMethod = fragmentMethod;
	}
	
	public JPeptideSpectrumMatch getPSM(){
		return this.jPeptideSpectrumMatch;
	}
	
	public void setPSM(JPeptideSpectrumMatch j){
		this.jPeptideSpectrumMatch = j;
	}

	/**
	 * Filters out the peaks with the lowest intensities until only 
	 * ticCutoffPercentage of the total ion current remains
	 * It's used by MVH
	 * @return The number of removed peaks
	 */
	public int filterByTIC(){
		int rpn=0;
		// Sort peak list in descending order based on intensity.
		// Use a multimap because multiple peaks can have the same intensity.
		// at least 7 fragment ion peaks.
		// Three class: 1, 2, 4 peaks
		if(this.jSpectrum.getPeaks().size()>7){
			double totalIonCurrent = this.jSpectrum.getTotalIonCurrent();
			double relativeIntensity=0.0;
			this.jSpectrum.sortPeaksByIntensity();
			ArrayList< JPeak> tmpJPeaks = new ArrayList<JPeak>();
			// sort peaks by intensity: from low to high
			for(int i=this.jSpectrum.getPeaks().size()-1;i>=0;i--){
				relativeIntensity+=this.jSpectrum.getPeaks().get(i).getIntensity()/totalIonCurrent;
				if(relativeIntensity<=this.ticCutoffPercentage){
					tmpJPeaks.add(this.jSpectrum.getPeaks().get(i));
				}else{
					rpn++;
				}
			}
			this.jSpectrum.setPeaks(tmpJPeaks);
		}
		return rpn;
	}


	public void setMaxFragmentChargeState(int maxFragmentChargeState) {
		this.maxFragmentChargeState = maxFragmentChargeState;
	}


	public int filterPeakByRange(double lowBoundmz,double upBoundmz){
		int rpn = 0;
		Iterator<JPeak> pListIterator = this.jSpectrum.getPeaks().iterator();  
		while(pListIterator.hasNext()){
			JPeak jPeak = pListIterator.next();
			if(jPeak.getMz() < lowBoundmz || jPeak.getMz() > upBoundmz){
				rpn++;
				pListIterator.remove();
			}
		}
		return rpn;
	}


	

	public void getSpectrumAnnotation() throws IllegalArgumentException{

		PeptideSpectrumAnnotator peptideSpectrumAnnotator = new PeptideSpectrumAnnotator();

		//String spectrumKey = spectrum.getSpectrumTitle();
		int charge = spectrum.getPrecursor().possibleCharges[0];
		//PeptideAssumption peptideAssumption = new PeptideAssumption(objPeptide, 1, 1, charge, 1.0);
		PeptideAssumption peptideAssumption = new PeptideAssumption(objPeptide, charge);
		SpecificAnnotationParameters specificAnnotationPreferences = new SpecificAnnotationParameters();

		HashSet<Integer> charges = new HashSet<>(4);
		int precursorCharge = peptideAssumption.getIdentificationCharge();
		if (precursorCharge == 1) {
			charges.add(precursorCharge);
		} else {
			for (int c = 1; c < precursorCharge; c++) {
				charges.add(c);
			}
		}
		specificAnnotationPreferences.setSelectedCharges(charges);

		specificAnnotationPreferences.addIonType(Ion.IonType.PEPTIDE_FRAGMENT_ION, PeptideFragmentIon.B_ION);
		specificAnnotationPreferences.addIonType(Ion.IonType.PEPTIDE_FRAGMENT_ION, PeptideFragmentIon.Y_ION);
		specificAnnotationPreferences.setFragmentIonAccuracy(CParameter.itol);
		specificAnnotationPreferences.setFragmentIonPpm(false);
		specificAnnotationPreferences.setNeutralLossesAuto(false);
		specificAnnotationPreferences.clearNeutralLosses();
		// this is important
		specificAnnotationPreferences.setPrecursorCharge(precursorCharge);

		if(lossWaterNH3) {
			specificAnnotationPreferences.addNeutralLoss(NeutralLoss.H2O);
			specificAnnotationPreferences.addNeutralLoss(NeutralLoss.NH3);
		}

		AnnotationParameters annotationSettings = new AnnotationParameters();
		annotationSettings.setTiesResolution(SpectrumAnnotator.TiesResolution.mostIntense);
		annotationSettings.setFragmentIonAccuracy(CParameter.itol);
		annotationSettings.setFragmentIonPpm(false);
		annotationSettings.setIntensityLimit(intensityLimit);
		annotationSettings.setNeutralLossesSequenceAuto(false);

		annotationSettings.setIntensityThresholdType(AnnotationParameters.IntensityThresholdType.percentile);
		//annotationSettings.clearNeutralLosses();

		if(lossWaterNH3) {
			annotationSettings.addNeutralLoss(NeutralLoss.H2O);
			annotationSettings.addNeutralLoss(NeutralLoss.NH3);
		}


		ModificationParameters modificationParameters = new ModificationParameters();
		SequenceMatchingParameters sequenceMatchingParameters = new SequenceMatchingParameters();
		JSequenceProvider jSequenceProvider = new JSequenceProvider();

		IonMatch[] matches = peptideSpectrumAnnotator.getSpectrumAnnotation(annotationSettings,
				specificAnnotationPreferences,
				"",
				"",
				spectrum,
				objPeptide,
				modificationParameters,
				jSequenceProvider,
				sequenceMatchingParameters);


        ArrayList<IonMatch> ionMatches = new ArrayList<>();

		if ((matches==null || matches.length ==0) && PeptideSearchMT.debug) {
			System.err.println("No ions matched!");
		}

		if(matches != null && matches.length > 0){
			Collections.addAll(ionMatches, matches);
		}

        this.matches = ionMatches;
       
	}

	public Peptide getObjPeptide(){
		if(this.objPeptide==null){
			this.objPeptide = new Peptide( jPeptideSpectrumMatch.getPepSeq());
			ArrayList<ModificationMatch> modificationMatches = jPeptideSpectrumMatch.getModificationMatch();
			for(ModificationMatch match : modificationMatches){
				this.objPeptide.addVariableModification(match);
			}

		}
		return this.objPeptide;
	}


	public boolean isLossWaterNH3() {
		return lossWaterNH3;
	}

	public void setLossWaterNH3(boolean lossWaterNH3) {
		this.lossWaterNH3 = lossWaterNH3;
	}






}
