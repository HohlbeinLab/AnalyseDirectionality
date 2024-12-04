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
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.InvalidPathException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

import static com.wurgobes.RFT.AnalyzerFunctions.*;
import static com.wurgobes.RFT.util.SaveCSV;
import static com.wurgobes.RFT.util.addLutLegend;


@Plugin(type = Command.class, name = "Angle Analyzer", menuPath = "Plugins>Angle Analyzer>Analyze Angles", priority = Priority.HIGH)
public class RFT<T extends RealType<T>> implements Command {
    private static final Logger log = LoggerFactory.getLogger(RFT.class);

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
    Plot hist_alt = new Plot("Angles", "Angle", "Intensity (a.u)");
    ArrayList<Roi> rois = new ArrayList<>();
    Overlay overlay = new Overlay();

    int width;
    int height;
    int map_width = -1;

    ImagePlus mask = null;
    ImagePlus max_imp = null;

    Path save_path = null;

    RFTParameters params = new RFTParameters();
    RFTDialogue<T> rftDialogue = null;

    static String macro_params = null;

    @Override
    public void run() {
        logService.info("Loaded RFT V1.0");

        circularLut = new OwnColorTable(lutService, logService);
        circularLut.setLut("spectrum");

        orderLUT = new OwnColorTable(lutService, logService);
        orderLUT.setLut("winter");

        hist_alt.savePlotObjects();


        if (macro_params == null)
            macro_params = Macro.getOptions();


        if (macro_params != null && !macro_params.isEmpty()) {
            
            params.getMacroParameters(macro_params);
            rftDialogue = new RFTDialogue<>(this, params);
        } else {
            rftDialogue = new RFTDialogue<>(this, params);
            rftDialogue.showDialog();
        }
    }

    public void setupImage() {

        if (imp == null || !imp.isVisible())
            imp = WindowManager.getCurrentImage();



        logService.info("RFT: Processing Image: " + imp.getTitle());


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

    public int saveData(RFTParameters params) {
        if(params.save_string == null || params.save_string.isEmpty()) {
            if (params.macro_mode)
                logService.error("Please ensure that the save_path is set in the macro.");
            else
                logService.error("Please ensure that the save-path is set in the settings pane.");
        } else {
            try {
                if (!params.save_string.contains("\\input\\")){
                    params.save_string = Paths.get(params.save_string, "input").toAbsolutePath().toString();
                }

                Files.createDirectories(Paths.get(params.save_string));
                save_path = Paths.get(params.save_string, "window" + params.window + "_" + imp.getShortTitle() + ".csv");


                logService.info("Saving to " + save_path);
                SaveCSV(csv_data, new ArrayList<>(Arrays.asList("x", "y", "width", "height", "Max Index", "Mask Median", "Angle", "Relevance?", "Profile Data")), save_path);
                return 1;
            } catch (IOException | InvalidPathException e) {
                logService.error("Failed to save data to path: " + params.save_string);
                logService.error(e.getMessage());
            }
        }
        return 0;
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
        int result = 0;
        if(params.macro_mode || (params.scanning_range && params.scan_save) || params.runPython)
            result = saveData(params);

        if(params.runPython && result==1)
            python(params);


        logService.info("Angle Analyzing Done");
    }

    public void calculateAngles(RFTParameters params){
        csv_data.clear();
        angle_map.clear();
        Pair<Integer, ImageStack> pair = AnalyzerFunctions.run(imp, max_imp, mask, csv_data, angle_map, circularLut, params);
        map_width = pair.getLeft();
        fft_stack = pair.getRight();

        if(!params.macro_mode){


            addLutLegend(max_imp.getProcessor(), circularLut, "Angle", 1024, 0f, 180f);
            max_imp.repaintWindow();
        }

        if (params.vector_overlay)
            calculateVectorField(params);


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
        for(int i = params.start; i <= params.end; i+=params.step){
            // For each window collect the angle dist and the order graph
            params.window = i;
            runVector(params);
            if(!params.macro_mode) {
                if (angle_stack == null)
                    angle_stack = new ImageStack(hist_alt.getProcessor().getWidth(), hist_alt.getProcessor().getHeight());

                angle_stack.addSlice("window " + i, hist_alt.getProcessor().duplicate());
            }
        }

        if (!params.macro_mode){
            ImagePlus angle_window_sweep = new ImagePlus("Angle vs Window", angle_stack);
            angle_window_sweep.show();
        }

    }

    public void python(RFTParameters params){
        logService.info("Running Python");
        if (!params.pythonPath.contains("py")){
            logService.info("Could not find Python at " + params.pythonPath);
            return;
        }
        Path scriptPath;
        try {
            scriptPath = Paths.get(params.scriptPath, "gaussian_order.py");
        } catch (Exception e) {
            logService.info("Could not find Script at " + params.scriptPath);
            return;
        }
        if (! new File(scriptPath.toString()).exists()){
            logService.info("Could not find Script at " + scriptPath);
            return;
        }
        if (! new File(save_path.toString()).exists()){
            logService.info("Could not find csv at " + scriptPath);
            return;
        }
        try {
            ProcessBuilder pb;
            if (params.showGraphs)
                pb = new ProcessBuilder(params.pythonPath, scriptPath.toString(), "-show_graph" , "-absolute", "-c", params.save_string, params.python_arguments, "-f", "window" + params.window + "_" + imp.getShortTitle(), params.python_arguments);
            else
                pb = new ProcessBuilder(params.pythonPath, scriptPath.toString() , "-absolute", "-c", params.save_string, params.python_arguments, "-f", "window" + params.window + "_" + imp.getShortTitle());
            pb.redirectErrorStream(true);
            Process pc = pb.start();
            BufferedReader pythonInput = new BufferedReader(new InputStreamReader(pc.getInputStream()));

            while (pc.isAlive()){
                String s;
                while ((s = pythonInput.readLine()) != null) {
                    logService.info("[Python] " + s);
                }

            }

        }
        catch (IOException e) {
            System.out.println("Failed to execute python Script. Did you set py?");
            e.printStackTrace();
        }
    }

    public void recordMacroParameters() {
            rftDialogue.recordMacroParameters();
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



