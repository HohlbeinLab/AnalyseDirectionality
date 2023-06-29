package com.wurgobes.AngleAnalyzer;

/*
Angle Analyzer
(c) 2022 Martijn Gobes, Wageningen University.


This software is released under the GPL v3. You may copy, distribute and modify
the software as long as you track changes/dates in source files. Any
modifications to or software including (via compiler) GPL-licensed code
must also be made available under the GPL along with build & install instructions.
https://www.gnu.org/licenses/gpl-3.0.en.html

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */


import gui_orientation.Credits;
import gui_orientation.Help;
import gui_orientation.WalkBarOrientationJ;
import gui_orientation.components.GridPanel;
import gui_orientation.components.GridToolbar;
import gui_orientation.components.Settings;
import gui_orientation.components.SpinnerDouble;
import ij.*;
import ij.gui.*;
import ij.plugin.frame.Recorder;
import ij.process.ColorProcessor;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import net.imagej.ImageJ;
import net.imagej.lut.LUTService;
import net.imglib2.type.numeric.RealType;
import org.apache.commons.lang3.tuple.Pair;
import org.scijava.Priority;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;


import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

import static com.wurgobes.AngleAnalyzer.AnalyzerFunctions.*;
import static com.wurgobes.AngleAnalyzer.util.*;


@Plugin(type = Command.class, name = "Angle Analyzer", menuPath = "Plugins>Angle Analyzer>Analyze Angles", priority = Priority.HIGH)
public class AngleAnalyzer <T extends RealType<T>> implements Command {

    // The services are passed through from ImageJ automatically

    @Parameter
    private LogService logService;

    @Parameter
    private LUTService lutService;


    /** The ImagePlus this plugin operates on. */
    //@Parameter(label="Image to process")
    //protected Dataset dataset;

    //@Parameter(type = ItemIO.OUTPUT)
    //private Dataset result;

    private ImagePlus imp;

    /* private variables */
    private static boolean debugging = false;

    private OwnColorTable circularLut; // Class to load LUT's
    private OwnColorTable orderLUT; // Class to load LUT's

    /* Directionality stuff */

    //x, y, width, height, index, median, angle, FT data
    ArrayList<ArrayList<Double>> csv_data = new ArrayList<>();
    ArrayList<Boolean> sig_map = new ArrayList<>();
    ArrayList<Double> adjusted_stats = new ArrayList<>();
    ArrayList<Double> angle_map = new ArrayList<>();
    ImageStack fft_stack = null;
    Plot hist_alt = new Plot("Histogram", "Angle", "Value");
    Plot order_plot = new Plot("Order", "Neighbourhood Size", "Avg Order");
    ArrayList<Roi> rois = new ArrayList<>();
    Overlay overlay = new Overlay();

    ImageStack order_stack = null;
    ImagePlus order_imp = null;

    int width;
    int height;
    int map_width = -1;

    ImagePlus mask = null;
    ImagePlus max_imp = null;

    RAFTParameters params = new RAFTParameters();


    @Override
    public void run() {
        logService.info("Angle Analyzer 0.1");


        circularLut = new OwnColorTable(lutService);
        circularLut.setLut("physics.lut");

        orderLUT = new OwnColorTable(lutService);
        orderLUT.setLut("winter");

        if(imp == null)
            imp = WindowManager.getCurrentImage();



        width = imp.getWidth();
        height = imp.getHeight();

        hist_alt.savePlotObjects();
        order_plot.savePlotObjects();

        // Order Parameter
        order_stack = new ImageStack(width, height);


        // Mask
        mask = imp.duplicate();
        mask.setTitle("Mask");
        util.MakeMask(mask);
        ImageConverter converter = new ImageConverter(mask);
        converter.convertToGray8();

        imp.setOverlay(overlay);


        RAFTDialogue<T> raftDialogue = new RAFTDialogue<>(this, params);
        raftDialogue.showDialog();
        System.out.println("done");
	}

    public void saveData(){
        SaveCSV(csv_data, new ArrayList<>(Arrays.asList("x", "y", "width", "height", "Max Index", "Mask Median", "Angle", "Relevance?", "Profile Data")), Paths.get(".\\test.csv"));
    }

    public void runVector(RAFTParameters params){
        if(params.window % 2 == 0)
            params.window += 1;
        //Non stack version for now
        if(false && imp.hasImageStack()){
            ImageStack input = imp.getStack();
            ImageStack imageStack = new ImageStack();
            System.out.println(imageStack.getSize());
            for(int i=0; i < input.getSize(); i++){
                System.out.println("Processing slice " + i);
                IJ.showStatus(String.format("Processing %d/%d", i, input.getSize()));
                ImagePlus dummy = new ImagePlus("dummy", input.getProcessor(i+1));
                ImageProcessor allpeaks = util.obtainDistancesFFT(dummy, 0.005f);
                imageStack.addSlice(String.valueOf(i), allpeaks);
            }
            ImagePlus result = new ImagePlus("FFT", imageStack);
            result.show();
        } else {
            ColorProcessor max_ip = new ColorProcessor(width, height);
            if (max_imp == null)
                max_imp = new ImagePlus("Colored blocks", max_ip);
            else
                max_imp.setProcessor(max_ip);


            max_imp.show();

            calculateAngles(params);
        }
        params.firstResults = Boolean.TRUE;
    }

    public void calculateAngles(RAFTParameters params){
        Pair<Integer, ImageStack> pair = AnalyzerFunctions.run(imp, max_imp, mask, csv_data, angle_map, circularLut, params.overlap, params.buffer, params.window, height, width);
        map_width = pair.getLeft();
        fft_stack = pair.getRight();
        ImagePlus fft_imp = new ImagePlus("ffts", fft_stack);
        fft_imp.show();
        addLutLegend(max_imp.getProcessor(), circularLut, "Angle", 1024, 0f, 180f);
        max_imp.repaintWindow();

        calculateVectorField(params);
        calculateOrder(params);
        AngleGraph(params);
        hist_alt.show();
    }


    public void calculateVectorField(RAFTParameters params){

        adjusted_stats.clear();
        calc_adjusted_stats(adjusted_stats, csv_data);
        applyVectorField(params);
    }

    public void applyVectorField(RAFTParameters params){
        rois.clear();

        calcVectorMap(rois, csv_data, params.window, params.vector_length, params.vector_width, params.cutoff);
        applyOverlay(imp, overlay, rois);
    }

    public void calculateOrder(RAFTParameters params){
        sig_map.clear();
        order_plot.restorePlotObjects();
        calcSigMap(csv_data, sig_map, params.intensity_cutoff);
        order_stack = order_parameter(csv_data, angle_map, sig_map, orderLUT, order_stack, order_plot, map_width, width, height, params.window);

        if(order_imp == null)
            order_imp = new ImagePlus("order_parameter", order_stack);

        if(order_imp.isVisible())
            order_imp.updateAndRepaintWindow();
        else
            order_imp.show();
    }

    public void toggleOverlay(){
        toggle_overlay(imp, overlay);
    }

    public void AngleGraph(RAFTParameters params){
        hist_alt.restorePlotObjects();
        createAngleGraph(fft_stack, csv_data, params.cutoff, hist_alt);
    }

    // Only run from the IDE
    public static void main(final String... arguments) {
        debugging = true;
        ImageJ ij = new ImageJ();
        ij.ui().showUI();

        //Dataset input = (Dataset) ij.io().open("test_stack.tif");
        //ij.ui().show(input);


        //ImagePlus imp = new ij.io.Opener().openImage("W:\\Data\\Microscopy\\RCM\\Test Data\\SPC Horizontal\\MAX_middle to bottom- bottom 30% 405 5% 561_1_MMStack_Pos0.ome-1-1.tif");
       // ImagePlus imp = new ij.io.Opener().openImage("W:\\Data\\Processing Testing\\MRI\\anisotropic_cropped.tif");
        //ImagePlus imp = new ij.io.Opener().openImage("W:\\Data\\Processing Testing\\Test Images\\tif\\45degright.tif");
        //ImagePlus imp = new ij.io.Opener().openImage("test.tif");
        //ImagePlus imp = new ij.io.Opener().openImage("C:\\Users\\marti\\Downloads\\Pretty pictures-20230626T220144Z-001\\Pretty pictures\\crop.tif");
        //ImagePlus imp = new ij.io.Opener().openImage("C:\\Users\\marti\\Downloads\\Pretty pictures-20230626T220144Z-001\\Pretty pictures\\600dpi (1).tif");
        //ImagePlus imp = new ij.io.Opener().openImage("C:\\Users\\marti\\Downloads\\MRI images for NWO\\CA_14_expr2_slice_5.tif");
        ImagePlus imp = new ij.io.Opener().openImage("C:\\Users\\marti\\Downloads\\MRI images for NWO\\CA_9_expr2_slice_5.tif");
        //ImagePlus imp = new ij.io.Opener().openImage("W:\\Data\\Microscopy\\Airyscan\\2022\\12\\Dead Stop SPC June 2022 11 cm\\10_crop.tif");
        //ImagePlus imp = new ij.io.Opener().openImage("W:\\Data\\Microscopy\\Airyscan\\2022\\12\\Dead Stop SPC June 2022 24 cm\\23_crop.tif");
        //ImagePlus imp = new ij.io.Opener().openImage("W:\\Data\\Microscopy\\Airyscan\\TVP_1\\TVP.tif");
        imp.show();
        //args = "D:\\Data\\Microscopy\\2022\\07\\8%561_40ms_MP3_1_RhB100x\\height slice_1\\slice_36_crop.tif";
        // "C:\\Users\\gobes001\\LocalSoftware\\AnalyseDirectionality\\test.tif";
        ij.command().run(AngleAnalyzer.class, true);
    }
}

class RAFTDialogue<T extends RealType<T>> extends JDialog implements ActionListener, ChangeListener, WindowListener, Runnable {

    private final AngleAnalyzer<T> angleAnalyzer;

    private final Settings settings				= new Settings("RAFT", IJ.getDirectory("plugins") + "RAFT.txt");
    protected RAFTParameters params;
    private final SpinnerDouble spnWindow					= new SpinnerDouble(50, 3, 100000, 1);
    private final SpinnerDouble spnOverlap		= new SpinnerDouble(0.75, 0, 1, 0.01);
    private final SpinnerDouble spnBuffer		= new SpinnerDouble(0, 0, 1, 0.01);
    private final SpinnerDouble spnVectorFieldLength		= new SpinnerDouble(70.0, 0, 10000, 1);
    private final SpinnerDouble spnVectorFieldWidth		= new SpinnerDouble(1.0, 0.1, 10, 0.1);
    private final SpinnerDouble spnCutoff		= new SpinnerDouble(2.0, 0, 10, 0.1);
    private final SpinnerDouble spnIntensityCutoff		= new SpinnerDouble(0, 0, 1, 0.01);

    private final JCheckBox	showVectorFieldOverlay	    = new JCheckBox("Overlay", true);
    protected WalkBarOrientationJ walk                  = new WalkBarOrientationJ();

    protected JButton				bnRun					= new JButton("Run");
    protected JButton				bnSave					= new JButton("Save Data");

    private enum Job {NONE, RUN, VECTOR_FIELD, SAVE, CUTOFF, ORDER, VECTOR_SCALE}
    private Job job = Job.NONE;
    private Thread					thread				= null;

    RAFTDialogue(AngleAnalyzer<T> angleAnalyzer, RAFTParameters params) {
        this.angleAnalyzer = angleAnalyzer;
        this.params = params;
        setTitle("RAFT");
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        if (e.getSource() == showVectorFieldOverlay)
            start(Job.VECTOR_FIELD);
        else if (e.getSource() == bnRun)
            start(Job.RUN);
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
        if (job == Job.VECTOR_FIELD && params.firstResults)
            angleAnalyzer.toggleOverlay();
        if (job == Job.VECTOR_SCALE && params.firstResults)
            angleAnalyzer.applyVectorField(params);

        if (job == Job.SAVE && params.firstResults)
            angleAnalyzer.saveData();
        if (job == Job.CUTOFF && params.firstResults) {
            angleAnalyzer.AngleGraph(params);
            angleAnalyzer.applyVectorField(params);
        }
        if (job == Job.ORDER && params.firstResults)
            angleAnalyzer.calculateOrder(params);


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

        pnMain.place(3, 0, pnVectors);
        pnMain.place(4, 0, pnCutoff);

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
        settings.record("spnVectorFieldScale", spnVectorFieldLength, "0.7");
        settings.record("spnVectorFieldScale", spnVectorFieldWidth, "1.0");


        settings.loadRecordedItems();
        setParameters();


    }

    public void getParameters() {
        params.buffer = spnBuffer.get();
        params.window = (int) spnWindow.get();
        params.overlap = spnOverlap.get();
        params.intensity_cutoff = spnIntensityCutoff.get();
        params.cutoff = spnCutoff.get();
        params.showVectorOverlay = showVectorFieldOverlay.isSelected();
        params.vector_length = spnVectorFieldLength.get()/100.0;
        params.vector_width =  spnVectorFieldWidth.get();

    }

    public void setParameters() {
        spnBuffer.set(params.buffer);
        spnWindow.set(params.window);
        showVectorFieldOverlay.setSelected(params.showVectorOverlay);
        spnOverlap.set(params.overlap);
        spnIntensityCutoff.set(params.intensity_cutoff);
        spnCutoff.set(params.cutoff);
    }

    private void recordMacroParameters() {
        if (!Recorder.record)
            return;
        String options = "";
        String plugin = "RAFT";

        options += "cutoff=" + spnCutoff.get() + " ";
        options += "window=" + spnWindow.get() + " ";
        options += "intensity_cutoff=" + spnIntensityCutoff.get() + " ";
        options += "overlap=" + spnOverlap.get() + " ";
        options += "buffer=" + spnBuffer.get() + " ";

        options += "vectorlength=" + spnVectorFieldLength.get() + " ";
        options += "vectorwidth=" + spnVectorFieldWidth.get() + " ";
        options += params.showVectorOverlay ? "vectoroverlay=on " : "vectoroverlay=off ";

        Recorder.record("run", plugin, options);
    }

    public RAFTParameters getSettingParameters() {
        return params;
    }

}

class RAFTParameters {

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

    public boolean firstResults = Boolean.FALSE;

    public Color vector_color = new Color(255, 227, 0);



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
