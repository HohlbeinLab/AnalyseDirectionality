package com.wurgobes.RFT;

import com.wurgobes.RFT.gui_components.*;
import ij.IJ;
import ij.gui.GUI;
import ij.plugin.frame.Recorder;
import net.imglib2.type.numeric.RealType;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;

class RFTDialogue<T extends RealType<T>> extends JDialog implements ActionListener, ChangeListener, WindowListener, Runnable {

    private final RFT<T> RFT;

    private final Settings settings = new Settings("RFT", IJ.getDirectory("plugins") + "RFT.txt");
    protected final RFTParameters params;
    private final SpinnerInteger spnWindow = new SpinnerInteger(75, 3, 100000, 1);
    private final SpinnerDouble spnOverlap = new SpinnerDouble(0.75, 0, 1, 0.01);
    private final SpinnerInteger spnBuffer = new SpinnerInteger(0, 0, 100000, 1);
    private final SpinnerDouble spnVectorFieldLength = new SpinnerDouble(100.0, 1.0, 10000, 1.0);
    private final SpinnerDouble spnVectorFieldWidth = new SpinnerDouble(3.0, 1.0, 10, 1.0);
    private final SpinnerDouble spnCutoff = new SpinnerDouble(2.0, 0, 10, 0.1);
    private final SpinnerDouble spnIntensityCutoff = new SpinnerDouble(0, 0, 1, 0.01);
    private final SpinnerInteger spnScanStart = new SpinnerInteger(21, 3, 100000, 1);
    private final SpinnerInteger spnScanEnd = new SpinnerInteger(61, 3, 100000, 1);
    private final SpinnerInteger spnScanStep = new SpinnerInteger(2, 2, 100000, 1);

    private final JTextField txtSavePath = new JTextField("save");
    private final JTextField txtPython = new JTextField("py");
    private final JTextField txtCore = new JTextField("core");
    private final JTextField pythonParams = new JTextField("params");
    private final JButton savePathButton = new JButton("Choose...");
    private final JButton pythonPathButton = new JButton("Choose...");
    private final JButton corePathButton = new JButton("Choose...");

    private final JCheckBox showVectorFieldOverlay = new JCheckBox("Overlay", true);
    private final JCheckBox saveDuringScan = new JCheckBox("Save Scans", true);
    private final JCheckBox runPython = new JCheckBox("Run Python", true);
    private final JCheckBox showGraph = new JCheckBox("Show Graph", false);
    protected WalkBarOrientationJ walk = new WalkBarOrientationJ();

    protected JButton bnRun = new JButton("Run");
    protected JButton bnSave = new JButton("Save Data");
    protected JButton bnScan = new JButton("Window Size Scan");
    protected JButton bnPython = new JButton("Save & Run Python");


    private enum Job {NONE, RUN, VECTOR_FIELD, SAVE, CUTOFF, ORDER, VECTOR_SCALE, SCAN, PYTHON, SAVE_PATH, PYTHON_PATH, CORE_PATH}

    private Job job = Job.NONE;
    private Thread thread = null;

    JFileChooser fc = new JFileChooser();

    RFTDialogue(RFT<T> RFT, RFTParameters params) {
        this.RFT = RFT;
        this.params = params;
        fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        setTitle("RFT");
        if (params.macro_mode) {
            setParameters();
            if(params.scanning_range)
                start(Job.SCAN);
            else
                start(Job.RUN);
        }

    }

    @Override
    public void actionPerformed(ActionEvent e) {

        if (e.getSource() == showVectorFieldOverlay)
            start(Job.VECTOR_FIELD);
        else if (e.getSource() == bnRun)
            start(Job.RUN);
        else if (e.getSource() == bnScan)
            start(Job.SCAN);
        else if (e.getSource() == bnSave)
            start(Job.SAVE);
        else if (e.getSource() == bnPython)
            start(Job.PYTHON);
        else if (e.getSource() == savePathButton)
            start(Job.SAVE_PATH);
        else if (e.getSource() == pythonPathButton)
            start(Job.PYTHON_PATH);
        else if (e.getSource() == corePathButton)
            start(Job.CORE_PATH);
        else if (e.getSource() == txtCore)
            params.scriptPath = txtCore.getText();
        else if (e.getSource() == txtSavePath)
            params.save_string = txtSavePath.getText();
        else if (e.getSource() == txtPython)
            params.pythonPath = txtPython.getText();
    }

    private void start(Job job) {
        if (thread != null)
            return;
        this.job = job;
        thread = new Thread(this);
        thread.setPriority(Thread.MIN_PRIORITY);
        thread.start();
    }

    @Override
    public void windowOpened(WindowEvent e) {

    }

    @Override
    public void windowClosing(WindowEvent e) {
        settings.storeRecordedItems();
        dispose();
    }

    @Override
    public void windowClosed(WindowEvent e) {

    }

    @Override
    public void windowIconified(WindowEvent e) {

    }

    @Override
    public void windowDeiconified(WindowEvent e) {

    }

    @Override
    public void windowActivated(WindowEvent e) {

    }

    @Override
    public void windowDeactivated(WindowEvent e) {
    }

    @Override
    public void run() {
        getParameters();
        walk.reset();

        if (job == Job.RUN)
            RFT.runVector(params);
        if (job == Job.SCAN)
            RFT.scanWindow(params);
        if (job == Job.VECTOR_FIELD && params.firstResults)
            RFT.toggleOverlay();
        if (job == Job.VECTOR_SCALE && params.firstResults)
            RFT.applyVectorField(params);


        if (job == Job.SAVE && params.firstResults)
            RFT.saveData(params);

        if (job == Job.CUTOFF && params.firstResults) {
            RFT.AngleGraph(params);
            RFT.applyVectorField(params);
        }
        if (job == Job.ORDER && params.firstResults) {
            RFT.AngleGraph(params);
            RFT.applyVectorField(params);
        }
        if (job == Job.PYTHON) {
            RFT.saveData(params);
            RFT.python(params);
        }
        if (job == Job.SAVE_PATH || job == Job.PYTHON_PATH || job == Job.CORE_PATH){

            if (fc.showDialog(this, "Select Path...") == JFileChooser.APPROVE_OPTION)
                switch(job){
                    case SAVE_PATH:
                        txtSavePath.setText(fc.getSelectedFile().getAbsolutePath());
                        break;
                    case PYTHON_PATH:
                        txtPython.setText(fc.getSelectedFile().getAbsolutePath());
                        break;
                    case CORE_PATH:
                        txtCore.setText(fc.getSelectedFile().getAbsolutePath());
                        break;
                }
            getParameters();
        }


        walk.finish();

        thread = null;
    }

    @Override
    public void stateChanged(ChangeEvent e) {
        if (e.getSource() == spnCutoff)
            start(Job.CUTOFF);
        else if (e.getSource() == spnVectorFieldLength || e.getSource() == spnVectorFieldWidth)
            start(Job.VECTOR_SCALE);
        else if (e.getSource() == spnIntensityCutoff)
            start(Job.ORDER);
    }

    public void showDialog() {
        // Panel Tensor
        GridToolbar pnOptions = new GridToolbar(false, 2);
        pnOptions.place(0, 0, new JLabel("Processing Window"));
        pnOptions.place(0, 2, spnWindow);
        pnOptions.place(0, 3, new JLabel("pixels"));
        pnOptions.place(1, 0, new JLabel("Overlap"));
        pnOptions.place(1, 2, spnOverlap);
        pnOptions.place(2, 0, new JLabel("Buffer"));
        pnOptions.place(2, 2, spnBuffer);
        pnOptions.place(2, 3, new JLabel("pixels"));
        pnOptions.place(3, 2, runPython);


        GridPanel pnMain1 = new GridPanel("Processing", 2);
        pnMain1.place(0, 0, pnOptions);
        pnMain1.place(2, 0, bnRun);
        pnMain1.place(3, 0, bnSave);
        pnMain1.place(4, 0, bnPython);

        GridPanel pnMain = new GridPanel(false);
        pnMain.place(0, 0, pnMain1);

        GridPanel pnVectors = new GridPanel("Vector Field");
        pnVectors.place(2, 0, new JLabel("Vector Length"));
        pnVectors.place(2, 1, spnVectorFieldLength);
        pnVectors.place(2, 2, new JLabel("%"));
        pnVectors.place(3, 0, new JLabel("Vector Width"));
        pnVectors.place(3, 1, spnVectorFieldWidth);
        pnVectors.place(3, 2, new JLabel("pixels"));
        pnVectors.place(6, 1, showVectorFieldOverlay);

        GridPanel pnScan = new GridPanel("Scan");
        pnScan.place(2, 0, new JLabel("Window Start"));
        pnScan.place(2, 1, spnScanStart);
        pnScan.place(2, 2, new JLabel("pixels"));
        pnScan.place(3, 0, new JLabel("Window Step"));
        pnScan.place(3, 1, spnScanStep);
        pnScan.place(3, 2, new JLabel("pixels"));
        pnScan.place(4, 0, new JLabel("Window End"));
        pnScan.place(4, 1, spnScanEnd);
        pnScan.place(4, 2, new JLabel("pixels"));
        pnScan.place(5, 1, bnScan);
        pnScan.place(6, 1, saveDuringScan);

        GridPanel pnCutoff = new GridPanel("Cutoffs");
        pnCutoff.place(2, 0, new JLabel("std Cutoff"));
        pnCutoff.place(2, 1, spnCutoff);
        pnCutoff.place(3, 0, new JLabel("Intensity Cutoff"));
        pnCutoff.place(3, 1, spnIntensityCutoff);


        showVectorFieldOverlay.addActionListener(this);
        spnVectorFieldLength.addChangeListener(this);
        spnVectorFieldWidth.addChangeListener(this);
        spnCutoff.addChangeListener(this);
        spnIntensityCutoff.addChangeListener(this);
        spnIntensityCutoff.addChangeListener(this);
        runPython.addActionListener(this);
        txtPython.addActionListener(this);
        txtSavePath.addActionListener(this);
        txtCore.addActionListener(this);
        runPython.addActionListener(this);
        savePathButton.addActionListener(this);
        corePathButton.addActionListener(this);
        pythonPathButton.addActionListener(this);
        txtCore.addActionListener(this);
        txtSavePath.addActionListener(this);

        pnMain.place(3, 0, pnVectors);
        pnMain.place(4, 0, pnCutoff);
        pnMain.place(5, 0, pnScan);

        GridPanel pnMainSettings = new GridPanel(false);
        GridPanel pnSettings = new GridPanel("Settings", 2);
        GridPanel pnOtherOptions = new GridPanel("Settings");
        pnOtherOptions.place(2, 0, new JLabel("Save Path"));
        pnOtherOptions.place(2, 1, txtSavePath);
        pnOtherOptions.place(2, 2, savePathButton);

        GridPanel pnPythonOptions = new GridPanel("Python");
        pnPythonOptions.place(2, 0, new JLabel("Python Path"));
        pnPythonOptions.place(2, 1, txtPython);
        pnPythonOptions.place(2, 2, pythonPathButton);
        pnPythonOptions.place(3, 0, new JLabel("Script Path"));
        pnPythonOptions.place(3, 1, txtCore);
        pnPythonOptions.place(3, 2, corePathButton);
        pnPythonOptions.place(4, 0, showGraph);
        pnPythonOptions.place(5, 0, new JLabel("Python Parameters"));
        pnPythonOptions.place(5, 1, pythonParams);

        pnSettings.place(0, 0, pnOtherOptions);
        pnSettings.place(1, 0, pnPythonOptions);

        pnMainSettings.place(2, 0, pnSettings);

        Help help = new Help();
        help.setPreferredSize(pnMain.getSize());

        Credits credits = new Credits();
        credits.setPreferredSize(pnMain.getSize());

        pnMainSettings.setPreferredSize(pnMain.getSize());

        JPanel pnHelp = new JPanel();
        pnHelp.setLayout(new BoxLayout(pnHelp, BoxLayout.PAGE_AXIS));

        JPanel pnHelpAdvanced = new JPanel(new BorderLayout());
        pnHelpAdvanced.add(pnHelp, BorderLayout.SOUTH);
        pnHelpAdvanced.add(help.getPane(), BorderLayout.CENTER);

        GridPanel pn = new GridPanel(false, 4);
        JTabbedPane tab = new JTabbedPane();

        tab.add("Processing", pnMain);
        tab.add("Settings", pnMainSettings);
        tab.add("Help", pnHelpAdvanced);
        tab.add("Credits", credits.getPane());

        pn.place(0, 0, tab);
        pn.place(1, 0, walk);

        // Listener
        walk.getButtonClose().addActionListener(this);
        bnRun.addActionListener(this);
        bnSave.addActionListener(this);
        bnScan.addActionListener(this);
        bnPython.addActionListener(this);

        // Finalize
        addWindowListener(this);
        getContentPane().add(pn);
        pack();
        setResizable(false);
        GUI.center(this);
        setVisible(true);

        settings.record("spnWindow", spnWindow, "50");
        settings.record("spnBuffer", spnBuffer, "0");
        settings.record("spnOverlap", spnOverlap, "0.75");
        settings.record("spnCutoff", spnCutoff, "2");
        settings.record("spnIntensityCutoff", spnIntensityCutoff, "0.2");
        settings.record("spnVectorFieldLength", spnVectorFieldLength, "100");
        settings.record("spnVectorFieldScale", spnVectorFieldWidth, "3.0");
        settings.record("start", spnScanStart, "17");
        settings.record("end", spnScanEnd, "51");
        settings.record("step", spnScanStep, "2");
        settings.record("scan_save", saveDuringScan, false);

        settings.record("savePath", txtSavePath, "save");
        settings.record("PythonPath", txtPython, "py");
        settings.record("ScriptPath", txtCore, "script");
        settings.record("ShowGraphs", showGraph, false);
        settings.record("autoPython", runPython, false);
        settings.record("pythonParams", pythonParams, "");

        settings.loadRecordedItems();
        getParameters();
    }

    public void getParameters() {
        params.buffer = spnBuffer.get();
        params.window = spnWindow.get();
        params.overlap = spnOverlap.get();
        params.intensity_cutoff = spnIntensityCutoff.get();
        params.cutoff = spnCutoff.get();
        params.vector_overlay = showVectorFieldOverlay.isSelected();
        params.vector_length = spnVectorFieldLength.get() / 100.0;
        params.vector_width = spnVectorFieldWidth.get();
        params.start = spnScanStart.get();
        params.end = spnScanEnd.get();
        params.step = spnScanStep.get();
        params.scan_save = saveDuringScan.isSelected();

        params.save_string = txtSavePath.getText();
        params.runPython = runPython.isSelected();
        params.pythonPath = txtPython.getText();
        params.scriptPath = txtCore.getText();
        params.showGraphs = showGraph.isSelected();
        params.python_arguments = pythonParams.getText();
    }

    public void setParameters() {
        spnBuffer.set(params.buffer);
        spnWindow.set(params.window);
        showVectorFieldOverlay.setSelected(params.vector_overlay);
        spnOverlap.set(params.overlap);
        spnIntensityCutoff.set(params.intensity_cutoff);
        spnCutoff.set(params.cutoff);

        spnScanStart.set(params.start);
        spnScanEnd.set(params.end);
        spnScanStep.set(params.step);

        spnVectorFieldLength.set(params.vector_length*100);
        spnVectorFieldWidth.set(params.vector_width);
        saveDuringScan.setSelected(params.scan_save);

        runPython.setSelected(params.runPython);
        txtSavePath.setText(params.save_string);
        txtPython.setText(params.pythonPath);
        txtCore.setText(params.scriptPath);
        showGraph.setSelected(params.showGraphs);
        pythonParams.setText(params.python_arguments);
    }

    public void recordMacroParameters() {
        if (!Recorder.record)
            return;

        String options = "";
        String plugin = "RFT";

        options += "cutoff=" + spnCutoff.get() + " ";
        options += "window=" + spnWindow.get() + " ";
        options += "intensity_cutoff=" + spnIntensityCutoff.get() + " ";
        options += "overlap=" + spnOverlap.get() + " ";
        options += "buffer=" + spnBuffer.get() + " ";

        options += "save_path=" + txtSavePath.getText() + " ";

        options += "vector_length=" + spnVectorFieldLength.get() + " ";
        options += "vector_width=" + spnVectorFieldWidth.get() + " ";
        options += params.vector_overlay ? "vector_overlay=True " : "vector_overlay=False ";

        if (saveDuringScan.isSelected()){
            options += "scan_save=True ";
            options += "start=" + spnScanStart.get() + " ";
            options += "end=" + spnScanEnd.get() + " ";
            options += "step=" + spnScanStep.get() + " ";
        }

        if (runPython.isSelected()){
            options += "run_python=True ";
            options += "python_path=" + txtPython.getText() + " ";
            options += "script_path=" + txtCore.getText() + " ";
            options += "python_arguments=" + params.python_arguments + " ";
            options += "show_graphs=" + showGraph.isSelected() + " ";
        }


        Recorder.record("run", plugin, options);
    }

    public RFTParameters getSettingParameters() {
        return params;
    }

}

