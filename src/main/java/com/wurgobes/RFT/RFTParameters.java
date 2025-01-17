package com.wurgobes.RFT;

import ij.Macro;

import java.awt.*;

public class RFTParameters {

    public boolean vector_overlay = true;
    public boolean macro_mode = false;

    public int buffer				= 0;
    public int				window				= 50;
    public double				cutoff			= 2;
    public double				overlap             = 0.75;
    public double				intensity_cutoff             = 0.2;
    public double vector_length = 1; //length
    public double vector_width = 3 ; // pixel width

    public int width;
    public int height;

    public boolean firstResults = Boolean.FALSE;

    public Color vector_color = new Color(255, 227, 0);
    public int start = 50;
    public int end = 300;
    public int step = 50;

    public boolean scanning_range = Boolean.FALSE;
    public boolean scan_save = Boolean.TRUE;

    public String save_string = null;
    public String pythonPath = null;
    public String scriptPath = null;
    public boolean showGraphs = Boolean.FALSE;
    public boolean runPython = Boolean.FALSE;

    public String python_arguments = "";

    public void getMacroParameters(String options) {
        macro_mode = true;
        //Own stuff
        buffer = Integer.parseInt(Macro.getValue(options, "buffer", "0"));
        window = Integer.parseInt(Macro.getValue(options, "window", "50"));
        cutoff = Double.parseDouble(Macro.getValue(options, "cutoff", "2.0"));
        overlap = Double.parseDouble(Macro.getValue(options, "overlap", "0.75"));
        intensity_cutoff = Double.parseDouble(Macro.getValue(options, "intensity_cutoff", "0.2"));
        start = Integer.parseInt(Macro.getValue(options, "start", "50"));
        end = Integer.parseInt(Macro.getValue(options, "end", "300"));
        step = Integer.parseInt(Macro.getValue(options, "step", "50"));
        scanning_range = Boolean.parseBoolean(Macro.getValue(options, "scanning_range", "False"));
        scan_save = Boolean.parseBoolean(Macro.getValue(options, "scan_save", "True"));

        // Vector Field
        vector_overlay = Macro.getValue(options, "vector_overlay", "off").equals("on");
        vector_length = Double.parseDouble(Macro.getValue(options, "vector_length", "1"));
        vector_width = Integer.parseInt(Macro.getValue(options, "vector_width", "3"));

        // Saving
        save_string = Macro.getValue(options, "save_path", null);

        // Python
        pythonPath = Macro.getValue(options, "python_path", null);
        scriptPath = Macro.getValue(options, "script_path", null);
        python_arguments = Macro.getValue(options, "python_arguments", null);

        showGraphs = Boolean.parseBoolean(Macro.getValue(options, "show_graphs", "False"));
        runPython = Boolean.parseBoolean(Macro.getValue(options, "run_python", "False"));
    }



}
