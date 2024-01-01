package zhanglab.inputconvertor.env;

public class InputConvertorConstants {

	// ## Generated features from Input Convertor
	public static final String IC_RT_FIELD_NAME = "ic_observed_rt";
	public static final String IC_TITLE_FIELD_NAME = "ic_title";
	public static final String IC_SCAN_NUM_FIELD_NAME = "ic_scan_num";
	public static final String IC_CHARGE_FIELD_NAME = "ic_charge";
	public static final String IC_PEPTIDE_FIELD_NAME = "ic_peptide";
	public static final String IC_INFERREND_PEPTIDE_FIELD_NAME = "InferredPeptide";
	
	public static final String IC_SA_FIELD_NAME = "ic_SA";
	public static final String IC_DELTA_RT_FIELD_NAME = "ic_abs(deltaRT)";
	public static final String IC_PPM_FIELD_NAME = "ic_ppm";
	
	public static final String IC_READ_FIELD_NAME = "Reads";
	public static final String IC_MQ_FIELD_NAME = "MeanQScore";
	public static final String IC_DELTA_SCORE_FIELD_NAME = "DeltaScore";
	public static final String IC_SPEC_ID_FEILD_NAME = "SpecID";
	public static final String IC_GENOMIC_ID_FEILD_NAME = "GenomicID";
	public static final String IC_LABEL_FEILD_NAME = "Label";
	public static final String IC_IS_CANONICAL_FEILD_NAME = "isCanonical";
	
	
	// ## Tool sets
	public static final String CASANOVO = "Casanovo";
	public static final String PNOVO3	= "pNovo3";
	public static final String PEAKS	= "PEAKS";
	public static final String PXG		= "pXg";
	public static final String COMET	= "Comet";
	public static final String AUTORT	= "AutoRT";
	public static final String IC_CONVERTOR		= "ic_convertor";
	
	
	// ## AutoRT header
	public static final String AUTORT_HEADER_X	=	"x";
	public static final String AUTORT_HEADER_Y	=	"y";
	public static final int	   AUTORT_X_IDX		=	0;
	public static final int	   AUTORT_Y_IDX		=	1;
	public static final int	   AUTORT_PRED_Y_IDX	=	2;
	
	// ## MS2PIP header
	public static final String MS2PIP_HEADER_SPECID	=	"spec_id";
	public static final String MS2PIP_HEADER_MODIFICATIONS	=	"modifications";
	public static final String MS2PIP_HEADER_PEPTIDE	=	"peptide";
	public static final String MS2PIP_HEADER_CHARGE	=	"charge";
	
	// ##Comet
	public static final String COMET_PEPTIDE_FIELD_NAME	=	"modified_peptide";
	public static final String COMET_RT_FIELD_NAME	=	"retention_time_sec";
	
	
	public static int getFieldIndex (String[] fields, String tag) {
		int index = -1;
		
		for(int i=0; i<fields.length; i++) {
			if(fields[i].equalsIgnoreCase(tag)) {
				index = i;
				break;
			}
		}
		
		return index;
	}
}
