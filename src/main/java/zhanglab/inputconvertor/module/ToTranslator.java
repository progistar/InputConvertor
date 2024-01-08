package zhanglab.inputconvertor.module;

import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;

public interface ToTranslator {

	public void doTranslate (CommandLine cmd) throws IOException, ParseException;
}
