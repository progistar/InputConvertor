package zhanglab.inputconvertor.module;

import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;

public interface NetMHCpan {
	
	public void toNetMHCpanInputFormat (CommandLine cmd) throws IOException, ParseException;
	
	public void addNetMHCpanOutput (CommandLine cmd) throws IOException, ParseException;
}
