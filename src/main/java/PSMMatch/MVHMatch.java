package main.java.PSMMatch;

import com.compomics.util.experiment.biology.ions.Ion;
import com.compomics.util.experiment.biology.ions.IonFactory;
import com.compomics.util.experiment.biology.ions.NeutralLoss;
import com.compomics.util.experiment.biology.ions.impl.PeptideFragmentIon;
import com.compomics.util.experiment.identification.matches.IonMatch;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import com.compomics.util.parameters.identification.search.ModificationParameters;
import main.java.pg.CParameter;
import main.java.pg.JSequenceProvider;
import main.java.util.Cloger;
import org.apache.commons.math.util.MathUtils;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public class MVHMatch extends JMatch {
	

	private double mvh=0.0;
	private double mzFidelity=0.0;
	private int numIntensityClasses=3;
	private int numMzFidelityClasses=3;
	private double mzLowerBound=200.0;
	private double mzUpperBound=2000.0;
	
	private int matchedPeaks = 0;
	
	public MVHMatch(){
		super();
	}

	public MVHMatch(Spectrum spectrum, JSpectrum jSpectrum, JPeptideSpectrumMatch JPeptideSpectrumMatch){
		super(spectrum, jSpectrum, JPeptideSpectrumMatch);
		
	}
	

	public boolean classifyPeakIntensities(){
		boolean yes =false;
		// 1:2:4, 7 peaks
		int totalPeaks = super.getJSpectrum().getPeaks().size();
		if(totalPeaks>=7){
			yes = true;

			int npeak = (int) Math.floor(1.0*totalPeaks/7.0);

			super.getJSpectrum().sortPeaksByIntensity();
			int nclass1 = npeak;
			int nclass2 = 2*npeak;
			for(int i=totalPeaks-1;i>=totalPeaks-nclass1;i--){
				//System.out.println("i="+i+","+totalPeaks+","+nclass1);
				super.getJSpectrum().getPeaks().get(i).setIntenClass(0);
			}

			for(int i=totalPeaks-1-nclass1;i>=totalPeaks-nclass1-nclass2;i--){
				super.getJSpectrum().getPeaks().get(i).setIntenClass(1);
			}

			for(int i=totalPeaks-1-nclass1-nclass2;i>=0;i--){
				super.getJSpectrum().getPeaks().get(i).setIntenClass(2);
			}

		}
		return yes;
	}


	public boolean SpectrumPreProcessing() throws IllegalArgumentException, FileNotFoundException, ClassNotFoundException, IOException, InterruptedException, SQLException{
		boolean yes = true;
		super.filterPeakByRange(this.mzLowerBound,this.mzUpperBound);
		super.removeParentPeak();
		super.filterByTIC();
		super.removeIsotopes();
		super.filterByPeakCount(300);
		if(this.classifyPeakIntensities()){
			yes=true;
		}else{
			yes=false;
		}
		return yes;
	}
	
	public double lnCombin(int n, int k){
		if(n<0 || k <0){
			return 0;
		}

		if(n<k){
			System.err.println("n should be greater than k: n="+n+",k="+k);
		}

		return MathUtils.factorialLog(n) - MathUtils.factorialLog(n-k) - MathUtils.factorialLog(k);   
	}
	

	public ArrayList<Double> removeFuzzyValues(HashSet<Double> va, double tol){

		Object[] val = va.toArray();
		Arrays.sort(val);
		ArrayList<Double> tmpVa = new ArrayList<Double>();
		double a1 = (Double) val[0];
		for(int i=1;i<val.length;i++){
			double a2 = (Double) val[i];
			if(a2 - a1 >= tol){
				tmpVa.add(a1);
			}
			a1 = a2;
		}
		tmpVa.add(a1);
		return tmpVa;
	}

	public void calcMVH() throws IllegalArgumentException, ClassNotFoundException, IOException, InterruptedException, SQLException{

		super.getJSpectrum().resetPeaks();
		if(!this.SpectrumPreProcessing()){
			return;
		}
		this.mvh=0.0;
		this.mzFidelity=0.0;

		int mvhKey[] =new int[this.numIntensityClasses+1];
		int mzFidelityKey[] = new int[this.numMzFidelityClasses+1];
		double[] mzFidelityThresholds = new double[this.numMzFidelityClasses];
		double minMzFidelityClassCount = (Math.pow(2.0, this.numMzFidelityClasses)-1)/(2-1);
		double lastI = 0;
		for(int i=0;i<this.numMzFidelityClasses-1;i++){
			mzFidelityThresholds[i] = CParameter.itol*(Math.pow(2, i)+lastI)/minMzFidelityClassCount;
			lastI=mzFidelityThresholds[i];
		}
		mzFidelityThresholds[this.numMzFidelityClasses-1] = CParameter.itol;
		if(JPeptideSpectrumMatch.outputDetail){
			for(int i=0;i<mzFidelityThresholds.length;i++){
				Cloger.getInstance().logger.info("mzFidenlity class "+i+":"+mzFidelityThresholds[i]);
			}
		}

		// new function
		ModificationParameters modificationParameters = new ModificationParameters();
		SequenceMatchingParameters sequenceMatchingParameters = new SequenceMatchingParameters();
		JSequenceProvider sequenceProvider = new JSequenceProvider();
		HashMap<Integer, HashMap<Integer, ArrayList<Ion>>> seqIonsMap = IonFactory.getInstance().getFragmentIons(super.getObjPeptide(),modificationParameters,sequenceProvider,sequenceMatchingParameters);
		ArrayList<Ion> seqIons = new ArrayList<Ion>();
		for(int iion : seqIonsMap.keySet()){
			HashMap<Integer, ArrayList<Ion>> tmpIons = seqIonsMap.get(iion);
			for(int jion: tmpIons.keySet()){
				seqIons.addAll(tmpIons.get(jion));
			}
			
		}

		int totalPredictedPeaks = 0;
		HashSet<Double> mzPredictUniqueSet = new HashSet<Double>();
		if(JPeptideSpectrumMatch.outputDetail){
			Cloger.getInstance().logger.info("");
		}
		ArrayList<Integer> chargeCandidate = new ArrayList<>();
		int maxConsiderFragCharge = 2;
		if(this.jPeptideSpectrumMatch.getCharge()==2){
			maxConsiderFragCharge = 1;
		}
		for(int ch=1;ch<=maxConsiderFragCharge;ch++){
			chargeCandidate.add(ch);
		}
		int jj=0;
		for(int i=0;i<seqIons.size();i++){
			Ion ion = seqIons.get(i);

			if(ion.getType()!= Ion.IonType.PEPTIDE_FRAGMENT_ION){
				continue;
			}

			PeptideFragmentIon peptideFragmentIon = (PeptideFragmentIon) ion;

			if(ion.getNeutralLosses().length>=2){
				continue;
			}

			if(!isLossWaterNH3() && ion.getNeutralLosses().length>=1){

				boolean findLoss = false;
				for(NeutralLoss nLoss:ion.getNeutralLosses()){
					if(nLoss.isSameAs(NeutralLoss.H2O) || nLoss.isSameAs(NeutralLoss.NH3)){
						findLoss = true;
						break;
					}
				}
				if(findLoss){
					continue;
				}
			}

			if(super.getFragmentMethod().equalsIgnoreCase("cid")){
				if(ion.getSubType()!=PeptideFragmentIon.B_ION && ion.getSubType()!=PeptideFragmentIon.Y_ION ){

					continue;
				}
				
			}else{
				System.err.println("Currently, we only support cid!");
				System.exit(0);
			}

			for(int j:chargeCandidate){

				double mz = ion.getTheoreticMz(j);
				if(mz >= this.mzLowerBound && mz <= this.mzUpperBound){
					mzPredictUniqueSet.add(mz);
					jj++;

					if(JPeptideSpectrumMatch.outputDetail){

						Cloger.getInstance().logger.info("jj="+jj+"\t"+ion.getName()+"\t"+peptideFragmentIon.getNumber()+"\t"+mz+"\t"+ion.getSubTypeAsString()+"\t"+ion.getNeutralLossesAsString()+"\t"+ion.getTypeAsString());
					}
				}
			}
		}

		totalPredictedPeaks = removeFuzzyValues(mzPredictUniqueSet, CParameter.itol).size();

		if(JPeptideSpectrumMatch.outputDetail){
			Cloger.getInstance().logger.info("totalPredictedPeaks="+totalPredictedPeaks);
		}

		HashMap< Double,JPeak> jPeaksSet = new HashMap<Double, JPeak>();

		int intenClassCounts[] = new int[this.numIntensityClasses+1];

		for(JPeak jPeak:super.getJSpectrum().getPeaks()){
			jPeaksSet.put(jPeak.getMz(), jPeak);

			intenClassCounts[jPeak.getIntenClass()]++;
		}
		double fragMassError = CParameter.itol;
		int totalPeakSpace = (int) ( this.mzUpperBound-this.mzLowerBound);
		intenClassCounts[this.numIntensityClasses] = (int) (Math.round(totalPeakSpace/(fragMassError*2.0)) - super.getJSpectrum().getPeaks().size());

		HashSet<Double> mzUniqueSet = new HashSet<Double>();
		for(IonMatch ionMatch: super.getIonMatches()){
			
			if(ionMatch.peakMz < this.mzLowerBound || ionMatch.peakMz > this.mzUpperBound){
				continue; 
			}

			if(ionMatch.ion.getSubType() != PeptideFragmentIon.B_ION && ionMatch.ion.getSubType() != PeptideFragmentIon.Y_ION){
				continue;
			}

			if(ionMatch.ion.hasNeutralLosses() && ionMatch.ion.getNeutralLosses().length>=2){
				continue;
			}
			
			double mz = ionMatch.peakMz;
			if(!jPeaksSet.containsKey(mz)){
				continue;
			}

			if(mzUniqueSet.contains(mz)){
				continue;
			}else{
				mzUniqueSet.add(mz);
			}

			mvhKey[jPeaksSet.get(mz).getIntenClass()]++;

			double mzError = Math.abs(ionMatch.getAbsoluteError());
			mzFidelityKey[classifyError(mzError, mzFidelityThresholds)]++;

		}

		matchedPeaks = mzUniqueSet.size();

		mvhKey[this.numIntensityClasses] = totalPredictedPeaks-matchedPeaks;
		mzFidelityKey[this.numMzFidelityClasses] = totalPredictedPeaks-matchedPeaks; ;

		int totalPeakBins = intenClassCounts[this.numIntensityClasses] + super.getJSpectrum().getPeaks().size();

		this.mvh=0.0;

		
		for(int i=0;i<=this.numIntensityClasses;i++){
			this.mvh += lnCombin(intenClassCounts[i],mvhKey[i]);
			if(JPeptideSpectrumMatch.outputDetail){
				Cloger.getInstance().logger.info("i="+i+","+intenClassCounts[i]+"\t"+mvhKey[i]);
			}
		}
		if(JPeptideSpectrumMatch.outputDetail){
			Cloger.getInstance().logger.info("total="+totalPeakBins+"\t"+totalPredictedPeaks);
		}
		this.mvh -= lnCombin( totalPeakBins, totalPredictedPeaks );
		this.mvh = 0.0 - this.mvh;


	}
	

	public int classifyError(double mzError, double[] mzFidelityThresholds){
		
		for(int i=0;i<mzFidelityThresholds.length;i++){
			if(Math.abs(mzError) <= mzFidelityThresholds[i]){
				return i;
			}
		}
		return mzFidelityThresholds.length -1;
	}

	public double getMvh() {
		return mvh;
	}

	public void setMvh(double mvh) {
		this.mvh = mvh;
	}

	public double getMzFidelity() {
		return mzFidelity;
	}

	public void setMzFidelity(double mzFidelity) {
		this.mzFidelity = mzFidelity;
	}

	public int getMatchedPeaks() {
		return matchedPeaks;
	}


}
