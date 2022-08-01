package main.java.PSMMatch;

import com.compomics.util.experiment.biology.ions.Ion;
import com.compomics.util.experiment.biology.ions.NeutralLoss;
import com.compomics.util.experiment.biology.ions.impl.PeptideFragmentIon;
import com.compomics.util.experiment.biology.modifications.Modification;
import com.compomics.util.experiment.biology.modifications.ModificationCategory;
import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import com.compomics.util.experiment.biology.modifications.ModificationType;
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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;


public final class JPeptideSpectrumMatch {

	public static boolean outputDetail = false;

	private String SpectrumID;

	private int SpectrumIndex;

	private String pepSeq;

	private ArrayList<JModification> modificationList = new ArrayList<JModification>();
	
	private int rank;
	private int charge;
	private double calculatedMassToCharge;
	private double experimentalMassToCharge;
	private double evalue;
	private ArrayList<String> proteins = new ArrayList<String>();
	private boolean isDecoy=false;
	private double rt=-1.0;


	private double annotationLevel = 0.02;


	public JPeptideSpectrumMatch(){
		
	}


	public ArrayList<ModificationMatch> getModificationMatch() {
		ArrayList<ModificationMatch> modificationMatches = new ArrayList<ModificationMatch>();
		ModificationFactory ptmFactory = ModificationFactory.getInstance();
		if (modificationList.size() >= 1) {
			for (int i = 0; i < modificationList.size(); i++) {
				JModification jModification = modificationList.get(i);
				double modMassDelta = Double.valueOf(String.format("%.4f", jModification.getModMassDelta()));
				String ptmName = "CMSPTM:" + jModification.getResidue() + ":" + modMassDelta;
				if (!ptmFactory.containsModification(ptmName)) {
					ArrayList<String> residuesArray = new ArrayList<String>();
					residuesArray.add(jModification.getResidue());
					Modification ptm = null;
					if (jModification.getResidue().equalsIgnoreCase("N-term")) {
						//ptm = new Modification(PTM.MODNP, ptmName, modMassDelta, new ArrayList<String>());
						ptm = new Modification(ModificationType.modn_peptide, ptmName, modMassDelta, null, ModificationCategory.Other);
					} else if (jModification.getResidue().equalsIgnoreCase("C-term")) {
						ptm = new Modification(ModificationType.modc_peptide, ptmName, modMassDelta, null, ModificationCategory.Other);
					} else {
						ptm = new Modification(ModificationType.modaa, ptmName, modMassDelta, residuesArray, ModificationCategory.Other);
					}
					ptmFactory.addUserModification(ptm);
				}
				ModificationMatch mm = null;
				if (jModification.getResidue().equalsIgnoreCase("N-term")) {
					mm = new ModificationMatch(ptmName, 1);
				} else if (jModification.getResidue().equalsIgnoreCase("C-term")) {
					mm = new ModificationMatch(ptmName, this.getPepSeq().length());
				} else {
					mm = new ModificationMatch(ptmName, jModification.getModLocation());
				}

				modificationMatches.add(mm);
			}

		}
		return (modificationMatches);
	}
	

	public double getMass(){
		Peptide peptideTmp  = new Peptide( getPepSeq());
		getModificationMatch();
		return peptideTmp.getMass();
	}

	
	public String getSpectrumID() {
		return SpectrumID;
	}


	public void setSpectrumID(String spectrumID) {
		SpectrumID = spectrumID;
	}


	public int getSpectrumIndex() {
		return SpectrumIndex;
	}



	public String getPepSeq() {
		return pepSeq;
	}


	public void setPepSeq(String pepSeq) {
		this.pepSeq = pepSeq;
	}



	public int getRank() {
		return rank;
	}


	public void setRank(int rank) {
		this.rank = rank;
	}


	public int getCharge() {
		return charge;
	}


	public void setCharge(int charge) {
		this.charge = charge;
	}


	public double getCalculatedMassToCharge() {
		return calculatedMassToCharge;
	}


	public void setCalculatedMassToCharge(double calculatedMassToCharge) {
		this.calculatedMassToCharge = calculatedMassToCharge;
	}


	public double getExperimentalMassToCharge() {
		return experimentalMassToCharge;
	}


	public void setExperimentalMassToCharge(double experimentalMassToCharge) {
		this.experimentalMassToCharge = experimentalMassToCharge;
	}


	public double getEvalue() {
		return evalue;
	}


	public void setEvalue(double evalue) {
		this.evalue = evalue;
	}

	public ArrayList<String> getProteins() {
		return proteins;
	}



	public boolean isDecoy() {
		return isDecoy;
	}

	public void setDecoy(boolean isDecoy) {
		this.isDecoy = isDecoy;
	}
	
	public int getLabel(){
		if(this.isDecoy){
			return -1;
		}else{
			return 1;
		}
	}


	public Peptide getPeptide(){

		Peptide objPeptide = new Peptide( this.getPepSeq());
		ArrayList<ModificationMatch> modificationMatches = this.getModificationMatch();
		for(ModificationMatch match : modificationMatches){
			objPeptide.addVariableModification(match);
		}
		return objPeptide;
	}
	

	public String getPeakAnnotation(Spectrum spectrum, JSpectrum jSpectrum){
		jSpectrum.resetPeaks();
		Peptide objPeptide = this.getPeptide();
		PeptideSpectrumAnnotator peptideSpectrumAnnotator = new PeptideSpectrumAnnotator();

		//String spectrumKey = spectrum.getSpectrumKey();
		int charge = getCharge();
		PeptideAssumption peptideAssumption = new PeptideAssumption(objPeptide, charge);
		SpecificAnnotationParameters specificAnnotationPreferences = new SpecificAnnotationParameters();
		
		HashSet<Integer> charges = new HashSet<Integer>(4);
		int precursorCharge = peptideAssumption.getIdentificationCharge();
		if (precursorCharge == 1) {
			charges.add(precursorCharge);
		} else {
		    for (int c = 1; c < precursorCharge; c++) {
		    	charges.add(c);
		    }
		}
		specificAnnotationPreferences.setSelectedCharges(charges);

		specificAnnotationPreferences.addIonType(Ion.IonType.PEPTIDE_FRAGMENT_ION,PeptideFragmentIon.B_ION);
		specificAnnotationPreferences.addIonType(Ion.IonType.PEPTIDE_FRAGMENT_ION,PeptideFragmentIon.Y_ION);
		specificAnnotationPreferences.setFragmentIonAccuracy(CParameter.itol);
		specificAnnotationPreferences.setFragmentIonPpm(false);
		specificAnnotationPreferences.setNeutralLossesAuto(false);
		specificAnnotationPreferences.clearNeutralLosses();
		specificAnnotationPreferences.addNeutralLoss(NeutralLoss.H2O);
		specificAnnotationPreferences.addNeutralLoss(NeutralLoss.NH3);

		AnnotationParameters annotationSettings = new AnnotationParameters();
		annotationSettings.setTiesResolution(SpectrumAnnotator.TiesResolution.mostIntense);
		annotationSettings.setFragmentIonAccuracy(CParameter.itol);
		annotationSettings.setFragmentIonPpm(false);
		annotationSettings.setIntensityLimit(this.annotationLevel);

		annotationSettings.setIntensityThresholdType(AnnotationParameters.IntensityThresholdType.snp);

		ModificationParameters modificationParameters = new ModificationParameters();
		SequenceMatchingParameters sequenceMatchingParameters = new SequenceMatchingParameters();
		JSequenceProvider jSequenceProvider = new JSequenceProvider();
		IonMatch matches[] = peptideSpectrumAnnotator.getSpectrumAnnotation(annotationSettings, specificAnnotationPreferences, "","",spectrum, objPeptide,
				modificationParameters,
				jSequenceProvider,
				sequenceMatchingParameters);

		StringBuilder mz_all = new StringBuilder();
		StringBuilder int_all = new StringBuilder();
		StringBuilder mz_match = new StringBuilder();
		StringBuilder int_match = new StringBuilder();
		StringBuilder m_label = new StringBuilder();

		HashMap<Double,Double> validPeakMz = new HashMap<Double, Double>();
		for(int i=0;i<jSpectrum.getPeaks().size();i++){
			JPeak jPeak = jSpectrum.getPeaks().get(i);
			validPeakMz.put(jPeak.getMz(),jPeak.getIntensity());
			if(i==0){
				mz_all.append(jPeak.getMz());
				int_all.append(jPeak.getIntensity());
			}else{
				mz_all.append(";"+jPeak.getMz());
				int_all.append(";"+jPeak.getIntensity());
			}
		}

		ArrayList<IonMatch> ionMatches = new ArrayList<>();
		for(IonMatch match: matches){
			ionMatches.add(match);
		}

		Iterator<IonMatch> pListIterator = ionMatches.iterator();
		HashSet<Double> ionHitSet = new HashSet<Double>();
		while (pListIterator.hasNext()) {
			IonMatch ionMatch = pListIterator.next();
			if (!validPeakMz.containsKey(ionMatch.peakMz) || ionHitSet.contains(ionMatch.peakMz) || ionMatch.ion.getNeutralLosses().length>=2) {
				pListIterator.remove();
				
			}else{
				ionHitSet.add(ionMatch.peakMz);
			}
		}

		
		
		for(int i=0;i<ionMatches.size();i++){
			IonMatch ionMatch = ionMatches.get(i);
			PeptideFragmentIon fragmentIon = ((PeptideFragmentIon) ionMatch.ion);
			String label = ionMatch.ion.getSubTypeAsString()+fragmentIon.getNumber()+ionMatch.ion.getNeutralLossesAsString()+ionMatch.charge;
			if(i==0){
				mz_match.append(ionMatch.peakMz);
				int_match.append(ionMatch.peakIntensity);
				m_label.append(label);
			}else{
				mz_match.append(";"+ionMatch.peakMz);
				int_match.append(";"+ionMatch.peakIntensity);
				m_label.append(";"+label);
			}
			
		}

		//IonFactory.getInstance().getDefaultNeutralLosses().clear();

		StringBuilder outline = new StringBuilder();
		outline.append(this.getSpectrumID()).append("\t");
		outline.append(this.getCalculatedMassToCharge()).append("\t");
		outline.append(this.getExperimentalMassToCharge()).append("\t");
		outline.append(this.getCharge()).append("\t");
		outline.append(this.getPepSeq()).append("\t");
		outline.append(m_label).append("\t");
		outline.append(mz_match).append("\t");
		outline.append(int_match).append("\t");
		outline.append(mz_all).append("\t");
		outline.append(int_all);
		
		return(outline.toString());
	}

	public double getRt() {
		return rt;
	}
	public void setRt(double rt) {
		this.rt = rt;
	}


}
