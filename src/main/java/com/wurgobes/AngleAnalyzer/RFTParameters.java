package com.wurgobes.AngleAnalyzer;

import ij.Macro;

import java.awt.*;

public class RFTParameters {

    public boolean				showVectorOverlay		= true;
    public boolean              macro_mode = false;
    public int					vectorGrid				= 10;
    public double				vectorScale				= 100;
    public int					vectorType				= 0;

    public int buffer				= 0;
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
    public int start = 50;
    public int end = 300;
    public int step = 50;

    public boolean scanning_range = Boolean.FALSE;

    public void getMacroParameters(String options) {
        macro_mode = true;
        //Own stuff
        buffer = Integer.parseInt(Macro.getValue(options, "buffer", "0"));
        window = Integer.parseInt(Macro.getValue(options, "window", "50"));
        cutoff = Integer.parseInt(Macro.getValue(options, "cutoff", "2"));
        overlap = Integer.parseInt(Macro.getValue(options, "overlap", "0.75"));
        intensity_cutoff = Integer.parseInt(Macro.getValue(options, "intensity_cutoff", "0.2"));
        start = Integer.parseInt(Macro.getValue(options, "start", "50"));
        end = Integer.parseInt(Macro.getValue(options, "end", "300"));
        step = Integer.parseInt(Macro.getValue(options, "step", "50"));
        scanning_range = Boolean.parseBoolean(Macro.getValue(options, "scanning_range", "False"));

        // Vector Field
        showVectorOverlay = Macro.getValue(options, "vectoroverlay", "on").equals("on");
        vectorGrid = Integer.parseInt(Macro.getValue(options, "vectorgrid", "10"));
        vectorScale = Double.parseDouble(Macro.getValue(options, "vectorscale", "100"));
        vectorType = Integer.parseInt(Macro.getValue(options, "vectortype", "0"));
    }



}
