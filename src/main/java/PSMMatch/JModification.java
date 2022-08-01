package main.java.PSMMatch;

/**
 * Modification Class
 *
 */
public class JModification {
	
	private int modLocation; //0-based
	private String residue;
	private double modMassDelta;
	public JModification(){
		
	}
	
	
	public int getModLocation() {
		return modLocation;
	}
	public void setModLocation(int modLocation) {
		this.modLocation = modLocation;
	}
	public String getResidue() {
		return residue;
	}
	public void setResidue(String residue) {
		this.residue = residue;
	}
	public double getModMassDelta() {
		return modMassDelta;
	}
	public void setModMassDelta(double modMassDelta) {
		this.modMassDelta = modMassDelta;
	}

}
