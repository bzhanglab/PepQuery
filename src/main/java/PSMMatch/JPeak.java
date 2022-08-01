package main.java.PSMMatch;


import com.compomics.util.experiment.biology.ions.impl.ElementaryIon;

/**
 * Fragment peak class
 *
 */
public class JPeak implements Cloneable{

	private double mz;



	/**
	 * Single charged mass
	 */
	private double mh;
	public double getMH() {
		if(mh > 0){
			return mh;
		}else{
			return getMH(this.mz,this.charge);
		}

	}

	public void setMH(double m) {
		this.mh = m;
	}

	public double getMH(double mz_value, int charge_value){
		double mass_value = mz_value * charge_value - charge_value * ElementaryIon.proton.getTheoreticMass() + ElementaryIon.proton.getTheoreticMass();
		if(mass_value < 0.0){
			mass_value = 0.0;
		}
		return mass_value;
	}




	private double intensity;

	public void setID(String ID) {
		this.ID = ID;
	}

	private String ID;

	public int getCharge() {
		return charge;
	}

	public void setCharge(int charge) {
		this.charge = charge;
	}

	private int charge = 0;

	/**
	 * Default is 0, real peak not complementary peak
	 */
	private int type = 0;

	public int getType() {
		return type;
	}

	public void setType(int type) {
		this.type = type;
	}


	// It's used for MHV scoring algorithm
	// The type of peak: 1, 2 or 3.
    // -1: no type
	private int intenClass=-1;

	public JPeak clone(){
		JPeak jPeak = null;
		try {
			jPeak = (JPeak) super.clone();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		return jPeak;
		
	}

	public JPeak(){

	}
	

	public JPeak(double mz,double intensity){
		this.mz = mz;
		this.intensity = intensity;
	}
	

	public JPeak(double mz,double intensity, int iclass){
		this.mz = mz;
		this.intensity = intensity;
		this.intenClass = iclass;
	}

    
    public JPeak(JPeak jPeak) {
        this.mz = jPeak.getMz();
        this.intensity = jPeak.getIntensity();
        this.intenClass = jPeak.getIntenClass();
    }

	public double getMz() {
		return mz;
	}

	public void setMz(double mz) {
		this.mz = mz;
	}

	public double getIntensity() {
		return intensity;
	}

	public void setIntensity(double intensity) {
		this.intensity = intensity;
	}


	public int getIntenClass() {
		return intenClass;
	}


	public void setIntenClass(int intenClass) {
		this.intenClass = intenClass;
	}


	public String getID(){
		return this.ID;
	}

	
}
