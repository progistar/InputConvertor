package zhanglab.inputconvertor.data;

public class ModificationTable {

	public static String COMET_M_OXI	= "M[15.9949]";
	public static String PEAKS_M_OXI	= "M(+15.99)";
	public static String PEAKS_C_CARBAM	= "C(+57.02)";
	
	
	public static String UNIMOD_M_OXI	= "M[UNIMOD:35]";
	public static String UNIMOD_C_CARBAM= "C[UNIMOD:4]";
	
	public static String PNOVO3_M_OXI	= "a";
	public static String PNOVO3_C_CARBAM	= "c";
	
	public static String MS2PIP_M_OXI = "Oxidation";
	public static String MS2PIP_C_CARBAM = "Carbamidomethyl";
	
	public static String IC_M_OXI	  = "M+15.995";
	public static String IC_C_CARBAM  = "C+57.021";

	public static String AUTORT_M_OXI = "1";
	public static String AUTORT_C_CARBAM = "2";
	
	
	public static final double OXIDATION_MASS = new ChemicalForm(0, 0, 1, 0, 0).getMass();
	public static final double CARBAM_MASS	  = new ChemicalForm(2, 3, 1, 1, 0).getMass();
}
