package zhanglab.inputconvertor.env;

public class ProteomeConstants {
	
	public static final double	Electron = 0.000549;
	public static final double	Hydrogen = 1.007825035;
	public static final double	Oxygen = 15.99491463;
	public static final double	Nitrogen = 14.00307401;
	public static final double 	Sulfur = 31.9720707;
	public static final double	Carbon = 12.0000000;
	public static final double	Phosphorus = 30.97376149;
	public static final double	Proton = 1.00727649;
	public static final double	HO = Hydrogen + Oxygen;	
	public static final double	H2O = Hydrogen*2 + Oxygen;	
	public static final double	NH3 = Hydrogen*3 + Nitrogen;		
	public static final double	IsotopeSpace = 1.00235;
	public static final double 	Sodium = 22.98976966;
	public static final double	Selenium = 79.9165213;
	public static final double	Ferrum = 55.9349418;
	public static final double	Potassium = 38.9637069;
	public static final double	Iodine = 126.90447;
	public static final double	B_ION_OFFSET = Proton;
	public static final double	Y_ION_OFFSET = H2O + Proton;
	
	public static final double	A_ION = Proton - (Carbon+Oxygen);
	public static final double	B_ION = Proton;
	public static final double	C_ION = Nitrogen + 3*Hydrogen + Proton;
	
	public static final double	X_ION = 2*Oxygen + Carbon + Proton;
	public static final double	Y_ION = H2O + Proton;
	public static final double	Z_ION = Oxygen - (Nitrogen+Hydrogen) + Proton;;
	
	public static final int UNIT_AS_DALTON = 0;
	public static final int UNIT_AS_PPM = 1;
	
	
	public static double ppmToleranceToDalton (double precursorMz, double ppmTolerance) {
		double dalton = ppmTolerance * precursorMz / Math.pow(10, 6);
		return dalton;
	}
}