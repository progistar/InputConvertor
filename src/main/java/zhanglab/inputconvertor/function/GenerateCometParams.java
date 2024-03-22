package zhanglab.inputconvertor.function;

import java.io.InputStream;

import org.apache.commons.cli.Options;

public class GenerateCometParams {

	//TODO: Generate comet params... one day!
	
	private static final String DEFAULT_DECOY_SEARCH = "1";
	private static final String DEFAULT_PEPTIDE_MASS_TOLERANCE= "20"; // ppm
	private static final String DEFAULT_ENZYME_NUMBER= "0"; // non-enzyme
	private static final String DEFAULT_ENZYME2_NUMBER= "0"; // non-enzyme
	private static final String DEFAULT_ENZYME_TERMINI= "2"; // fully, but it does not work for non-enzyme
	private static final String DEFAULT_ALLOWED_MISSED_CLEAVAGE= "2";
	
	private static final String DEFAULT_FRAGMENT_BIN_TOL= "0.02"; // for high res
	private static final String DEFAULT_FRAGMENT_BIN_OFFSET= "0.0"; // for high res
	private static final String DEFAULT_SAMPLE_ENZYME_NUMBER= "0"; // for non-enzyme
	
	private static final String DEFAULT_OUTPUT_SUFFIX= "MANDATORY"; // output suffix is a mandatory option.
	private static final String DEFAULT_ADD_NTERM_PEPTIDE= "0.0"; // None
	private static final String DEFAULT_ADD_C_CYSTEINIE= "0.0"; // None
	private static final String DEFAULT_ADD_K_LYSINE= "0.0"; // None

	public GenerateCometParams (String[] args) {
		
		Options options = new Options();
		
		// Options
		options.addOption("m", "module", true, "mode (ex> pxg_input, ms2pip_input, autort_input");
		
		options.addOption("d", "decoy", true, "0=no, 1=internal decoy concatenated (default), 2=internal decoy separate");
		options.addOption("t", "peptide_mass_tolerance", true, "peptide mass tolerance at ppm");
		options.addOption("E", "enzyme1", true, "");
		
		
		ClassLoader classLoader = GenerateCometParams.class.getClassLoader();
		System.out.println(classLoader.getName());
        // Read the resource using the ClassLoader
        try (InputStream inputStream = classLoader.getResourceAsStream("comet.params.new")) {
            if (inputStream != null) {
            	
            } else {
                System.err.println("Resource not found!");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
	}
	
}
