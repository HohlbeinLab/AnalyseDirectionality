package com.wurgobes.AngleAnalyzer;

import ij.gui.GenericDialog;


import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

@Plugin(type = Command.class, menuPath = "Plugins>Angle Analyzer>About Angle Analyzer...")
public class About implements Command {

    @Parameter
    private LogService log;

    @Override
    public void run() {
        log.info("Angle Analyzer About");
        String content = "Developed by Martijn Gobes at the Holhbein Lab.\nMore information can be found at https://github.com/HohlbeinLab/\n" +
                "Current Version: 0.1";

        GenericDialog gd = new GenericDialog("About");
        gd.addMessage(content);

        gd.showDialog();
    }


}
