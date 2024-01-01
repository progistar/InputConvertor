package zhanglab.inputconvertor.module;

import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;

public interface ToFeatures {
	
	public void toFeatures (CommandLine cmd) throws IOException, ParseException;
	
}
