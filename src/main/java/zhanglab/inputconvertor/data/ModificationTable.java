package zhanglab.inputconvertor.data;

public class ModificationTable {

	public static String COMET_M_OXI	= "M[15.9949]";
	public static String PEAKS_M_OXI	= "M(+15.99)";
	public static String PNOVO3_M_OXI	= "a";
	
	public static String MS2PIP_M_OXI = "Oxidation";
	public static String MS2PIP_C_CAM = "Carbamidomethyl";
	
	public static String IC_M_OXI	  = "M+15.995";

	public static String AUTORT_M_OXI = "1";
	
	
	public static final double OXIDATION_MASS = new ChemicalForm(0, 0, 1, 0, 0).getMass();
}
