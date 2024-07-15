package com.wurgobes.AngleAnalyzer;

import gui_orientation.Credits;
import gui_orientation.Help;
import gui_orientation.WalkBarOrientationJ;
import gui_orientation.components.*;
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

    private final AngleAnalyzer<T> angleAnalyzer;

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

    private final JCheckBox showVectorFieldOverlay = new JCheckBox("Overlay", true);
    private final JCheckBox saveDuringScan = new JCheckBox("Save Scans", true);
    protected WalkBarOrientationJ walk = new WalkBarOrientationJ();

    protected JButton bnRun = new JButton("Run");
    protected JButton bnSave = new JButton("Save Data");
    protected JButton bnScan = new JButton("Window Size Scan");


    private enum Job {NONE, RUN, VECTOR_FIELD, SAVE, CUTOFF, ORDER, VECTOR_SCALE, SCAN}

    private Job job = Job.NONE;
    private Thread thread = null;


    RFTDialogue(AngleAnalyzer<T> angleAnalyzer, RFTParameters params) {
        this.angleAnalyzer = angleAnalyzer;
        this.params = params;
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
            angleAnalyzer.runVector(params);
        if (job == Job.SCAN)
            angleAnalyzer.scanWindow(params);
        if (job == Job.VECTOR_FIELD && params.firstResults)
            angleAnalyzer.toggleOverlay();
        if (job == Job.VECTOR_SCALE && params.firstResults)
            angleAnalyzer.applyVectorField(params);


        if (job == Job.SAVE && params.firstResults)
            angleAnalyzer.saveData(params);

        if (job == Job.CUTOFF && params.firstResults) {
            angleAnalyzer.AngleGraph(params);
            angleAnalyzer.applyVectorField(params);
        }
        if (job == Job.ORDER && params.firstResults) {
            angleAnalyzer.AngleGraph(params);
            angleAnalyzer.applyVectorField(params);
            angleAnalyzer.calculateOrder(params);
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
        GridToolbar pnTensor = new GridToolbar(false, 2);
        pnTensor.place(0, 0, new JLabel("Processing Window"));
        pnTensor.place(0, 2, spnWindow);
        pnTensor.place(0, 3, new JLabel("pixel"));
        pnTensor.place(1, 0, new JLabel("Overlap"));
        pnTensor.place(1, 2, spnOverlap);
        pnTensor.place(2, 0, new JLabel("Buffer"));
        pnTensor.place(2, 2, spnBuffer);


        GridPanel pnMain1 = new GridPanel("Processing", 2);
        pnMain1.place(0, 0, pnTensor);
        //pnMain1.place(1, 0, pnFeatures);
        pnMain1.place(2, 0, bnRun);
        pnMain1.place(3, 0, bnSave);

        GridPanel pnMain = new GridPanel(false);
        pnMain.place(0, 0, pnMain1);

        GridPanel pnVectors = new GridPanel("Vector Field");
        pnVectors.place(2, 0, new JLabel("Vector Length(%)"));
        pnVectors.place(2, 1, spnVectorFieldLength);
        pnVectors.place(3, 0, new JLabel("Vector Width (px)"));
        pnVectors.place(3, 1, spnVectorFieldWidth);
        pnVectors.place(6, 1, showVectorFieldOverlay);

        GridPanel pnScan = new GridPanel("Scan");
        pnScan.place(2, 0, new JLabel("Window Start"));
        pnScan.place(2, 1, spnScanStart);
        pnScan.place(3, 0, new JLabel("Window Step"));
        pnScan.place(3, 1, spnScanStep);
        pnScan.place(4, 0, new JLabel("Window End"));
        pnScan.place(4, 1, spnScanEnd);
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

        pnMain.place(3, 0, pnVectors);
        pnMain.place(4, 0, pnCutoff);
        pnMain.place(5, 0, pnScan);

        Help help = new Help();
        help.setPreferredSize(pnMain.getSize());

        JPanel pnHelp = new JPanel();
        pnHelp.setLayout(new BoxLayout(pnHelp, BoxLayout.PAGE_AXIS));

        JPanel pnHelpAdvanced = new JPanel(new BorderLayout());
        pnHelpAdvanced.add(pnHelp, BorderLayout.SOUTH);
        pnHelpAdvanced.add(help.getPane(), BorderLayout.CENTER);

        GridPanel pn = new GridPanel(false, 4);
        JTabbedPane tab = new JTabbedPane();
        Credits credits = new Credits();
        credits.setPreferredSize(pnMain.getSize());
        tab.add("Processing", pnMain);
        tab.add("Help", pnHelpAdvanced);
        tab.add("Credits", credits.getPane());
        pn.place(0, 0, tab);
        pn.place(1, 0, walk);

        // Listener
        walk.getButtonClose().addActionListener(this);
        bnRun.addActionListener(this);
        bnSave.addActionListener(this);
        bnScan.addActionListener(this);

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
        settings.record("spnVectorFieldScale", spnVectorFieldLength, "100");
        settings.record("spnVectorFieldScale", spnVectorFieldWidth, "3.0");
        settings.record("start", spnScanStart, "17");
        settings.record("end", spnScanEnd, "51");
        settings.record("step", spnScanStep, "2");


        settings.loadRecordedItems();
        setParameters();


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
    }

    private void recordMacroParameters() {
        if (!Recorder.record)
            return;

        String options = "";
        String plugin = "RFT";

        options += "cutoff=" + spnCutoff.get() + " ";
        options += "window=" + spnWindow.get() + " ";
        options += "intensity_cutoff=" + spnIntensityCutoff.get() + " ";
        options += "overlap=" + spnOverlap.get() + " ";
        options += "buffer=" + spnBuffer.get() + " ";

        options += "vectorlength=" + spnVectorFieldLength.get() + " ";
        options += "vectorwidth=" + spnVectorFieldWidth.get() + " ";
        options += params.vector_overlay ? "vectoroverlay=True " : "vectoroverlay=False ";

        Recorder.record("run", plugin, options);
    }

    public RFTParameters getSettingParameters() {
        return params;
    }

}

