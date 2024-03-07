package com.wurgobes.AngleAnalyzer;

import ij.Macro;

import java.awt.*;

public class RFTParameters {

    public boolean vector_overlay = true;
    public boolean macro_mode = false;
    public boolean order = false;

    public int buffer				= 0;
    public int				window				= 50;
    public double				cutoff			= 2;
    public double				overlap             = 0.75;
    public double				intensity_cutoff             = 0.2;
    public double vector_length = 1; //length
    public int vector_width = 3 ; // pixel width

    public int width;
    public int height;

    public boolean firstResults = Boolean.FALSE;

    public Color vector_color = new Color(255, 227, 0);
    public int start = 50;
    public int end = 300;
    public int step = 50;

    public boolean scanning_range = Boolean.FALSE;

    public String save_string = null;

    public void getMacroParameters(String options) {
        macro_mode = true;
        //Own stuff
        buffer = Integer.parseInt(Macro.getValue(options, "buffer", "0"));
        window = Integer.parseInt(Macro.getValue(options, "window", "50"));
        cutoff = Integer.parseInt(Macro.getValue(options, "cutoff", "2"));
        overlap = Double.parseDouble(Macro.getValue(options, "overlap", "0.75"));
        intensity_cutoff = Double.parseDouble(Macro.getValue(options, "intensity_cutoff", "0.2"));
        start = Integer.parseInt(Macro.getValue(options, "start", "50"));
        end = Integer.parseInt(Macro.getValue(options, "end", "300"));
        step = Integer.parseInt(Macro.getValue(options, "step", "50"));
        scanning_range = Boolean.parseBoolean(Macro.getValue(options, "scanning_range", "False"));

        // Vector Field
        vector_overlay = Macro.getValue(options, "vectoroverlay", "off").equals("on");
        vector_length = Double.parseDouble(Macro.getValue(options, "vector_length", "1"));
        vector_width = Integer.parseInt(Macro.getValue(options, "vector_width", "3"));

        // Saving
        save_string = Macro.getValue(options, "path", null);

        order = Boolean.parseBoolean(Macro.getValue(options, "path", "false"));
    }



}
