package main.java.PSMMatch;


import com.compomics.util.experiment.biology.ions.impl.PeptideFragmentIon;
import com.compomics.util.experiment.identification.matches.IonMatch;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;

import main.java.pg.PeptideSearchMT;
import main.java.util.Cloger;
import org.apache.commons.math.util.MathUtils;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Scoring algorithm: Hyperscore
 *
 */
public final class HyperscoreMatch extends JMatch {

	/**
	 * The number of matched b ions
	 */
	private int count_b_ion = 0;
	/**
	 * The number of matched y ions
	 */
	private int count_y_ion = 0;
	private double hyperScore=0.0;
	private double dScale = 4.0;

	/**
	 * In order to improve speed, only need to generate one time.
	 */
	private static HashMap<Integer,Double> cacheFactorialValues = new HashMap<>();


	public HyperscoreMatch(Spectrum spectrum, JSpectrum jSpectrum, JPeptideSpectrumMatch JPeptideSpectrumMatch){
		super(spectrum, jSpectrum, JPeptideSpectrumMatch);
		
	}

	/**
	 * Pre-processing a spectrum
	 * @throws IllegalArgumentException
	 */
	public void SpectrumPreProcessing() throws IllegalArgumentException{
		//super.getJSpectrum().getPeaks().size()

		super.removeIsotopes();
		//remove ions near the parent ion m/z
		super.removeParentPeak();
		//remove low mass immonium ions
		super.removeLowMasses();
		super.removeLowIntensityPeak();
		super.cleanIsotopes();
		super.filterByPeakCount();
		super.normalizePeaks(100.0);
	}

	public void calcHyperScore() throws IllegalArgumentException{
		getJSpectrum().resetPeaks();
		this.hyperScore = 0.0;
		this.count_b_ion=0;
		this.count_y_ion=0;
		double dScore=0.0;
		// Pre-processing spectrum
		this.SpectrumPreProcessing();

		// jPeaksSet saves the valid peaks (not filtered by pre-processing step)
		HashMap< Double,Double> jPeaksSet = new HashMap<Double, Double>();
		super.getJSpectrum().sortPeaksByMZ();
		for(JPeak jPeak:super.getJSpectrum().getPeaks()){
			jPeaksSet.put(jPeak.getMz(), jPeak.getIntensity());
			if(JPeptideSpectrumMatch.outputDetail){
				Cloger.getInstance().logger.info("detect peak: "+jPeak.getMz()+"\t"+ jPeak.getIntensity());
			}
		}

		if(super.getFragmentMethod().equalsIgnoreCase("cid") ||
				super.getFragmentMethod().equalsIgnoreCase("hcd")){

			ArrayList<Double> bArrayList = new ArrayList<>();

			ArrayList<Double> yArrayList = new ArrayList<>();
			// It's possible that the same fragment ion has multiple matches (different ion types)
			HashSet<Double> mzUniqueSet = new HashSet<>();

			for(IonMatch ionMatch:super.getIonMatches()){
				PeptideFragmentIon peptideFragmentIon = (PeptideFragmentIon) ionMatch.ion;

				//if(peptideFragmentIon.getNeutralLossesAsString().contentEquals("-H2O-NH3")){
				if(ionMatch.ion.hasNeutralLosses() && ionMatch.ion.getNeutralLosses().length>=2){
					continue;
				}

				if(JPeptideSpectrumMatch.outputDetail){
					System.err.println("Loss: "+peptideFragmentIon.getNeutralLossesAsString()+", Type:"+ionMatch.getPeakAnnotation());
				}

				// If the peak is not filtered by the pre-processing step, then this peak will be used for scoring
				if(jPeaksSet.containsKey(ionMatch.peakMz)){

					if(ionMatch.ion.getSubType()==PeptideFragmentIon.B_ION){

						if(mzUniqueSet.contains(ionMatch.peakMz)){
							continue;
						}else{
							mzUniqueSet.add(ionMatch.peakMz);
						}
						bArrayList.add(jPeaksSet.get(ionMatch.peakMz));
						count_b_ion++;
						dScore+=jPeaksSet.get(ionMatch.peakMz);

						if(JPeptideSpectrumMatch.outputDetail){
							Cloger.getInstance().logger.info("b:"+super.getPSM().getSpectrumIndex()+"\t"+ionMatch.charge+"\t"+ionMatch.peakMz+"\t"+ionMatch.peakIntensity+"\t"+ionMatch.getPeakAnnotation());
						}

					}else if(ionMatch.ion.getSubType()==PeptideFragmentIon.Y_ION){
						if(mzUniqueSet.contains(ionMatch.peakMz)){
							continue;
						}else{
							mzUniqueSet.add(ionMatch.peakMz);
						}
						yArrayList.add(jPeaksSet.get(ionMatch.peakMz));
						count_y_ion++;
						dScore+=jPeaksSet.get(ionMatch.peakMz);

						if(JPeptideSpectrumMatch.outputDetail){
							Cloger.getInstance().logger.info("y:"+super.getPSM().getSpectrumIndex()+"\t"+ionMatch.charge+"\t"+ionMatch.peakMz+"\t"+ionMatch.peakIntensity+"\t"+ionMatch.getPeakAnnotation());
						}
					}else{

						if(JPeptideSpectrumMatch.outputDetail){
							Cloger.getInstance().logger.info("Matched ion:"+ionMatch.peakMz+"=>"+ionMatch.peakIntensity+", ("+ionMatch.getPeakAnnotation()+") has not been used!");
							Cloger.getInstance().logger.info("n:"+super.getPSM().getSpectrumIndex()+"\t"+ionMatch.peakMz+"\t"+ionMatch.peakIntensity+"\t"+ionMatch.getPeakAnnotation());
						}
					}
				}else{
					if(JPeptideSpectrumMatch.outputDetail){

						Cloger.getInstance().logger.info("e:"+super.getPSM().getSpectrumIndex()+"\t"+ionMatch.charge+"\t"+ionMatch.peakMz+"\t"+ionMatch.peakIntensity+"\t"+ionMatch.getPeakAnnotation());
					}
					
				}
				
			}


			if(this.count_b_ion>64){
				this.count_b_ion=64;
			}
			if(this.count_y_ion>64){
				this.count_y_ion=64;
			}

			if(JPeptideSpectrumMatch.outputDetail){
				Cloger.getInstance().logger.info("b:"+this.count_b_ion+"; "+"y: "+this.count_y_ion);
			}
			
			if(this.count_b_ion > 0 && this.count_y_ion > 0){
				this.hyperScore = Math.log10( dScore )+cacheFactorialValues.get(this.count_b_ion)+cacheFactorialValues.get(this.count_y_ion);
			}else if(this.count_b_ion==0 && this.count_y_ion >0 ){
				// no b ion
				this.hyperScore = Math.log10( dScore )+cacheFactorialValues.get(this.count_y_ion);
			}else if(this.count_y_ion==0 && this.count_b_ion >0 ){
				// no y ion
				this.hyperScore = Math.log10( dScore )+cacheFactorialValues.get(this.count_b_ion);
			}else{

				if(PeptideSearchMT.debug) {
					Cloger.getInstance().logger.debug("Don't find valid matched peak for spectrum:" + super.getMSnSpectrum().getSpectrumTitle());
				}

				this.hyperScore=0.0;
				return;
			}
			
			this.hyperScore = this.dScale*this.hyperScore;
			
		}else{

			Cloger.getInstance().logger.error("Currently, we only support cid/hcd!");
			System.exit(0);
		}
		
	}


	public double getHyperScore() {
		return hyperScore;
	}

	public void setHyperScore(double hyperScore) {
		this.hyperScore = hyperScore;
	}

	public int getCount_b_ion() {
		return count_b_ion;
	}


	public int getCount_y_ion() {
		return count_y_ion;
	}



	/**
	 * This is used to improve speed.
	 * @param maxValue The max value
	 */
	public static void generateFactorialValues(int maxValue){
		for(int i=1;i<=maxValue;i++){
			cacheFactorialValues.put(i,Math.log10( 1.0*MathUtils.factorialDouble(i)));
		}
		cacheFactorialValues.put(0,0.0);
	}
	
}
