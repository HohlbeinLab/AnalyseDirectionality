package com.wurgobes.RFT.gui_components;

/**
 * This class and all other classes in this folder were adapted from (<a href="https://github.com/Biomedical-Imaging-Group/OrientationJ">OrientationJ</a>)
 * Please refer to their GitHub page for further documentation.
 */
public class WalkBarOrientationJ extends WalkBar {

    /**
     * Constructor
     */
    public WalkBarOrientationJ() {
        super(Constants.copyright, true, false, true);
        fillAbout(
                Constants.softname + " " + Constants.version,
                Constants.date,
                "",
                Constants.author,
                "Biophysics, Wageningen University & Research, Wageningen, Netherlands",
                "",
                "Citation Here");

    }
}
