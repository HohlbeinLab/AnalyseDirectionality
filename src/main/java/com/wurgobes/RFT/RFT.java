package com.wurgobes.RFT;

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


import ij.*;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.Roi;
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

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

import static com.wurgobes.RFT.AnalyzerFunctions.*;
import static com.wurgobes.RFT.util.SaveCSV;
import static com.wurgobes.RFT.util.addLutLegend;


@Plugin(type = Command.class, name = "Angle Analyzer", menuPath = "Plugins>Angle Analyzer>Analyze Angles", priority = Priority.HIGH)
public class RFT<T extends RealType<T>> implements Command {

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
    ImagePlus fft_imp = null;
    Plot hist_alt = new Plot("Angles", "Angle", "Intensity (a.u)");
    Plot order_plot = new Plot("Order", "Neighbourhood Size", "Avg Order");
    ArrayList<Roi> rois = new ArrayList<>();
    Overlay overlay = new Overlay();


    ImagePlus order_imp = null;

    int width;
    int height;
    int map_width = -1;

    ImagePlus mask = null;
    ImagePlus sig_mask = null;
    ImagePlus max_imp = null;

    RFTParameters params = new RFTParameters();

    static String macro_params = null;

    @Override
    public void run() {
        logService.info("Loaded RFT V1.0");

        circularLut = new OwnColorTable(lutService, logService);
        circularLut.setLut("spectrum");

        orderLUT = new OwnColorTable(lutService, logService);
        orderLUT.setLut("winter");

        hist_alt.savePlotObjects();
        order_plot.savePlotObjects();


        if (macro_params == null)
            macro_params = Macro.getOptions();


        if (macro_params != null && !macro_params.isEmpty()) {
            params.getMacroParameters(macro_params);
            new RFTDialogue<>(this, params);
        } else {
            RFTDialogue<T> rftDialogue = new RFTDialogue<>(this, params);
            rftDialogue.showDialog();
        }
    }

    public void setupImage() {

        imp = WindowManager.getCurrentImage();

        logService.info("Angle Analyzer 1.0. Processing Image: " + imp.getTitle());


        width = imp.getWidth();
        height = imp.getHeight();

        params.width = width;
        params.height = height;

        // Mask
        mask = imp.duplicate();
        mask.setTitle("Mask");
        util.MakeMask(mask);
        ImageConverter converter = new ImageConverter(mask);
        converter.convertToGray8();

        imp.setOverlay(overlay);
    }

    public void saveData(RFTParameters params) {
        try {
            Path path;
            if(params.save_string == null) {
                path = Paths.get("C:\\Users\\gobes001\\LocalSoftware\\AnalyseDirectionality\\Python Scripts\\input\\window" + params.window + "_" + imp.getShortTitle() + ".csv");
            }
            else {
                Files.createDirectories(Paths.get(params.save_string));
                path = Paths.get(params.save_string, "window" + params.window + "_" + imp.getShortTitle() + ".csv");
            }

            logService.info("Saving to " + path);
            SaveCSV(csv_data, new ArrayList<>(Arrays.asList("x", "y", "width", "height", "Max Index", "Mask Median", "Angle", "Relevance?", "Profile Data")), path);
        } catch (IOException e) {
            logService.error("Failed to save data to path: " + params.save_string);
            throw new RuntimeException(e);
        }
    }

    public void runVector(RFTParameters params){
        if(params.window % 2 == 0)
            params.window += 1;

        setupImage();
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

            if(!params.macro_mode){
                ColorProcessor max_ip = new ColorProcessor(width, height);
                if (max_imp == null)
                    max_imp = new ImagePlus("Colored blocks", max_ip);
                else
                    max_imp.setProcessor(max_ip);


                max_imp.show();
            }

            calculateAngles(params);
        }
        params.firstResults = Boolean.TRUE;

        if(params.macro_mode || (params.scanning_range && params.scan_save)){
            saveData(params);

        }


        logService.info("Angle Analyzing Done");
    }

    public void calculateAngles(RFTParameters params){
        csv_data.clear();
        angle_map.clear();
        Pair<Integer, ImageStack> pair = AnalyzerFunctions.run(imp, max_imp, mask, csv_data, angle_map, circularLut, params);
        map_width = pair.getLeft();
        fft_stack = pair.getRight();

        if(!params.macro_mode){
            if(fft_imp == null)
                fft_imp = new ImagePlus("ffts", fft_stack);
            else
                fft_imp.setStack(fft_stack);

            addLutLegend(max_imp.getProcessor(), circularLut, "Angle", 1024, 0f, 180f);
            max_imp.repaintWindow();
        }

        if (params.vector_overlay)
            calculateVectorField(params);
        if (params.order)
            calculateOrder(params);

        if(!params.macro_mode) {
            AngleGraph(params);
            hist_alt.show();
        }
    }

    public void calculateVectorField(RFTParameters params){

        adjusted_stats.clear();
        calc_adjusted_stats(adjusted_stats, csv_data);
        applyVectorField(params);
    }

    public void applyVectorField(RFTParameters params){
        rois.clear();

        calcVectorMap(rois, csv_data, params);
        applyOverlay(imp, overlay, rois, params);
    }

    public void calculateOrder(RFTParameters params){
        sig_map.clear();
        order_plot.restorePlotObjects();
        ImageProcessor sig_ip = calcSigMap(csv_data, sig_map, params);
        ImageStack order_stack = order_parameter(csv_data, angle_map, sig_map, orderLUT, order_plot, map_width, width, height, params.window);


        if(!params.macro_mode){
            if(order_imp == null)
                order_imp = new ImagePlus("order_parameter", order_stack);
            else
                order_imp.setStack(order_stack);

            if(sig_mask == null)
                sig_mask = new ImagePlus("Intensity Cutoff Mask", sig_ip);
            else
                sig_mask.setProcessor(sig_ip);

            order_plot.setLimitsToFit(Boolean.TRUE);
            order_plot.show();

            if(order_imp.isVisible())
                order_imp.updateAndRepaintWindow();
            else
                order_imp.show();

            if(sig_mask.isVisible())
                sig_mask.updateAndRepaintWindow();
            else
                sig_mask.show();
        }
    }

    public void toggleOverlay(){
        toggle_overlay(imp, overlay);
    }

    public void AngleGraph(RFTParameters params){
        hist_alt.restorePlotObjects();
        createAngleGraph(fft_stack, csv_data, params.cutoff, hist_alt);
    }

    public void scanWindow(RFTParameters params){
        params.scanning_range = Boolean.TRUE;
        ImageStack angle_stack = null;
        ImageStack order_graph_stack = null;
        for(int i = params.start; i <= params.end; i+=params.step){
            // For each window collect the angle dist and the order graph
            params.window = i;
            runVector(params);
            if(!params.macro_mode) {
                if (angle_stack == null)
                    angle_stack = new ImageStack(hist_alt.getProcessor().getWidth(), hist_alt.getProcessor().getHeight());
                if (order_graph_stack == null)
                    order_graph_stack = new ImageStack(order_plot.getProcessor().getWidth(), order_plot.getProcessor().getHeight());

                angle_stack.addSlice("window " + i, hist_alt.getProcessor().duplicate());
                order_graph_stack.addSlice("window " + i, order_plot.getProcessor().duplicate());
            }
        }

        if (!params.macro_mode){
            ImagePlus angle_window_sweep = new ImagePlus("Angle vs Window", angle_stack);
            ImagePlus order_window_sweep = new ImagePlus("Order vs Window", order_graph_stack);
            angle_window_sweep.show();
            order_window_sweep.show();
        }

    }

    // Only run from the IDE
    public static void main(final String... arguments) {
        debugging = true;
        ImageJ ij = new ImageJ();
        ij.ui().showUI();

        //macro_params = "buffer=0 window=250 overlap=0.5";
        ImagePlus imp = new ij.io.Opener().openImage("C:\\Users\\gobes001\\source\\repos\\Scratch\\2D FFT images\\mosaic2.tif");


        imp.show();
        ij.command().run(RFT.class, true);
    }
}



