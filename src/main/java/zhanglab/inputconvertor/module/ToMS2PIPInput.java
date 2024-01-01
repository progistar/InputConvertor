package zhanglab.inputconvertor.module;

import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;

public interface ToMS2PIPInput {

	
	public void toMS2PIPInputFormat (CommandLine cmd) throws IOException, ParseException;
}
