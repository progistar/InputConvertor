package zhanglab.inputconvertor.env;

public class InputConvertorConstants {

	public static final String VERSION = "Version: v1.0.0a";
	
	// ## Generated features from Input Convertor
	public static final String IC_RT_FIELD_NAME = "ic_observed_rt";
	public static final String IC_TITLE_FIELD_NAME = "ic_title";
	public static final String IC_SCAN_NUM_FIELD_NAME = "ic_scan_num";
	public static final String IC_CHARGE_FIELD_NAME = "ic_charge";
	public static final String IC_PEPTIDE_FIELD_NAME = "ic_peptide";
	public static final String IC_SEARCH_SCORE_FIELD_NAME = "ic_search_score";
	public static final String IC_INFERRED_PEPTIDE_FIELD_NAME = "InferredPeptide";
	
	public static final String IC_SA_FIELD_NAME = "ic_SA";
	public static final String IC_DELTA_RT_FIELD_NAME = "ic_abs(deltaRT)";
	public static final String IC_PPM_FIELD_NAME = "ic_ppm";
	
	public static final String IC_READ_FIELD_NAME = "Reads";
	public static final String IC_MQ_FIELD_NAME = "MeanQScore";
	public static final String IC_DELTA_SCORE_FIELD_NAME = "DeltaScore";
	public static final String IC_SPEC_ID_FEILD_NAME = "SpecID";
	public static final String IC_GENOMIC_ID_FEILD_NAME = "GenomicID";
	public static final String IC_GENE_NAME_FEILD_NAME = "GeneName";
	public static final String IC_LABEL_FEILD_NAME = "Label";
	public static final String IC_IS_CANONICAL_FEILD_NAME = "IsCanonical";
	
	
	// ## Tool sets
	public static final String CASANOVO = "Casanovo";
	public static final String PNOVO3	= "pNovo3";
	public static final String PEAKS	= "PEAKS";
	public static final String PXG		= "pXg";
	public static final String COMET	= "Comet";
	public static final String AUTORT	= "AutoRT";
	public static final String MS2PIP	= "ms2pip";
	public static final String PROSIT	= "prosit";
	public static final String IC_CONVERTOR		= "ic_convertor";
	public static final String ARRIBA	= "Arriba";
	public static final String EXPERIMENTAL_SPECTRUM = "experimental_spectrum";
	
	
	// ## Casanovo header
	public static final String CASANOVO_PEPTIDE_FIELD_NAME = "sequence";
	public static final String CASANOVO_SPECTRA_REF_FIELD_NAME = "spectra_ref";
	public static final String CASANOVO_SCORE_FIELD_NAME = "search_engine_score[1]";
	
	// ## PEAKS header
	public static final String PEAKS_SOURCE_FILE_FIELD_NAME = "Source File";
	public static final String PEAKS_SCAN_FIELD_NAME = "Scan";
	public static final String PEAKS_PEPTIDE_FIELD_NAME = "Peptide";
	public static final String PEAKS_CHARGE_FIELD_NAME = "z";
	public static final String PEAKS_SCORE_FIELD_NAME = "ALC (%)";
	
	
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
	
	// ## Comet header
	public static final String COMET_PEPTIDE_FIELD_NAME	=	"modified_peptide";
	public static final String COMET_RT_FIELD_NAME	=	"retention_time_sec";
	public static final String COMET_MODIFICATION_FIELD_NAME = "modifications";
	public static final String COMET_PROTEIN_FIELD_NAME = "protein";
	public static final String COMET_SCAN_FIELD_NAME = "scan";
	public static final String COMET_CHARGE_FIELD_NAME = "charge";
	public static final String COMET_XCORR_FIELD_NAME = "xcorr";
	
	public static final String COMET_PIN_SPECID_FIELD_NAME = "specid";
	public static final String COMET_PIN_PEPLEN_FIELD_NAME = "peplen";
	public static final String COMET_PIN_PEPTIDE_FIELD_NAME = "peptide";
	public static final String COMET_PIN_PROTEIN_FIELD_NAME = "proteins";
	public static final String COMET_PIN_SCAN_FIELD_NAME = "scannr";
	
	public static final String COMET_PIN_DM_FIELD_NAME = "dm";
	public static final String COMET_PIN_ABSDM_FIELD_NAME = "absdm";
	public static final String COMET_PIN_EXPMASS_FIELD_NAME = "expmass";
	public static final String COMET_PIN_CALCMASS_FIELD_NAME = "calcmass";
	public static final String COMET_PIN_MASS_FIELD_NAME = "mass";
	public static final String COMET_PIN_ENZN_FIELD_NAME = "enzn";
	public static final String COMET_PIN_ENZC_FIELD_NAME = "enzc";
	public static final String COMET_PIN_ENZ_INT_FIELD_NAME = "enzint";
	
	// ## VEP header
	public static final String VEP_LOC_FIELD_NAME = "Location";
	public static final String VEP_REF_FIELD_NAME = "REF_ALLELE";
	public static final String VEP_ALT_FIELD_NAME = "Allele";
	public static final byte				WILD   = 0;
	public static final byte				SNP   = 1;
	public static final byte				MNP   = 2;
	public static final byte				INS   = 3;
	public static final byte				DEL   = 4;
	public static final byte				ALL_MUT   = 5;
	public static final String DELETION_MARK	  = "-";
	public static final String MUTATION_HEADER_ID = "@VR";
	
	// ## Arriba
	public static final String ARRIBA_PEPTIDE_FIELD_NAME	=	"peptide_sequence";
	public static final String ARRIBA_GENE1_FIELD_NAME		=	"#gene1";
	public static final String ARRIBA_GENE2_FIELD_NAME		=	"gene2";
	public static final String ARRIBA_BREAK_POINT1_FIELD_NAME		=	"breakpoint1";
	public static final String ARRIBA_BREAK_POINT2_FIELD_NAME		=	"breakpoint2";
	public static final String ARRIBA_HEADER_ID				=	"@FG";
	public static final String ARRIBA_FRAME_FIELD_NAME		=	"reading_frame";
	public static final String ARRIBA_STRAND1_FIELD_NAME	=	"strand1(gene/fusion)";
	public static final String ARRIBA_STRAND2_FIELD_NAME	=	"strand2(gene/fusion)";
	public static final String ARRIBA_TRANSCRIPT1_FIELD_NAME=	"transcript_id1";
	public static final String ARRIBA_TRANSCRIPT2_FIELD_NAME=	"transcript_id2";
	
	// ## StringTie
	public static final String STRINGTIE_HEADER_ID			=	"@ST";
	
	// ## IRFinder
	public static final String IRFINDER_HEADER_ID			=	"@IR";
	public static final String IRFINDER_CHR_FIELD_NAME		=	"Chr";
	public static final String IRFINDER_START_FIELD_NAME	=	"Start";
	public static final String IRFINDER_END_FIELD_NAME		=	"End";
	public static final String IRFINDER_NAME_FIELD_NAME		=	"Name";
	public static final String IRFINDER_STRAND_FIELD_NAME	=	"Strand";
	public static final String IRFINDER_WARNINGS_FIELD_NAME	=	"Warnings";
	
	// ## CIRIQuant
	public static final String CIRIQUANT_HEADER_ID			=	"@CR";
	
	// ## GTF Threeframe translation
	public static final String EXON_TRANSLATION_HEADER_ID	=	"@EX";
	
	// ## Reference
	public static final String REF_HEADER_ID				=	"@RF";
	public static final String NON_REF_HEADER_ID			=	"@NF";
	
	// ## Translation
	public static final int MIN_PEPT_LEN					=	7;
	
	
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
