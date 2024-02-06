package com.wurgobes.AngleAnalyzer;

import ij.Macro;

import java.awt.*;

public class RAFTParameters {

    public boolean				showVectorOverlay		= true;

    public int					vectorGrid				= 10;
    public double				vectorScale				= 100;
    public int					vectorType				= 0;

    public double				buffer				= 0.0;
    public int				window				= 50;
    public double				cutoff			= 2;
    public double				overlap             = 0.75;
    public double				intensity_cutoff             = 0.2;
    public double vector_length;
    public double vector_width;

    public int width;
    public int height;

    public boolean firstResults = Boolean.FALSE;

    public Color vector_color = new Color(255, 227, 0);
    public int start = 17;
    public int end = 101;


    public void getMacroParameters(String options) {
        //Own stuff
        buffer = Integer.parseInt(Macro.getValue(options, "buffer", "0"));
        window = Integer.parseInt(Macro.getValue(options, "window", "50"));
        cutoff = Integer.parseInt(Macro.getValue(options, "cutoff", "2"));
        overlap = Integer.parseInt(Macro.getValue(options, "overlap", "0.75"));
        intensity_cutoff = Integer.parseInt(Macro.getValue(options, "intensity_cutoff", "0.2"));

        // Vector Field
        showVectorOverlay = Macro.getValue(options, "vectoroverlay", "on").equals("on");
        vectorGrid = Integer.parseInt(Macro.getValue(options, "vectorgrid", "10"));
        vectorScale = Double.parseDouble(Macro.getValue(options, "vectorscale", "100"));
        vectorType = Integer.parseInt(Macro.getValue(options, "vectortype", "0"));
    }



}
