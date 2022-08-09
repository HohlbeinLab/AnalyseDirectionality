package com.wurgobes.AngleAnalyzer;

import gui_orientation.AnalysisDialog;
import gui_orientation.WalkBarOrientationJ;
import ij.plugin.PlugIn;
import orientation.*;
import orientation.imageware.ImageWare;

public class OrientationJ_Vector_Field implements PlugIn {

    public void run(String arg) {

        if (arg == null) {
            AnalysisDialog orientation = new AnalysisDialog(OrientationService.VECTORFIELD);
            orientation.showDialog();
        } else {
            OrientationParameters params = new OrientationParameters(OrientationService.VECTORFIELD);
            params.getMacroParameters(arg);
            ImageWare source = GroupImage.getCurrentImage();
            if (source == null) {
                return;
            }
            WalkBarOrientationJ walk = new WalkBarOrientationJ();
            OrientationProcess process = new OrientationProcess(walk, source, params);
            process.run();
            // Version 2.0.5: ansewer to pull request
            OrientationResults.show(process.getGroupImage(), params, 1);
        }


    }
}
