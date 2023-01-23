package com.wurgobes.AngleAnalyzer;

/*
Angle Analyzer
(c) 2022 Martijn Gobes, Wageningen University.


The distance and angle is also averaged.

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


import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.ProfilePlot;
import ij.process.ImageProcessor;
import net.imagej.ImageJ;
import net.imagej.lut.LUTService;
import net.imglib2.type.numeric.RealType;
import org.scijava.Priority;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import java.awt.*;
import java.io.IOException;
import java.util.Arrays;


@Plugin(type = Command.class, name = "Angle Analyzer", menuPath = "Plugins>Angle Analyzer>Analyze Angles", priority = Priority.HIGH)
public class AngleAnalyzer <T extends RealType<T>> implements Command {

    // The services are passed through from ImageJ automatically

    @Parameter
    private LogService logService;

    @Parameter
    private LUTService lutService;

    @Parameter(label = "Size of processing window")
    private float window = 150; //px

    /** The ImagePlus this plugin operates on. */
    //@Parameter(label="Image to process")
    //protected Dataset dataset;

    //@Parameter(type = ItemIO.OUTPUT)
    //private Dataset result;

    private ImagePlus imp;

    /* private variables */
    private static boolean debugging = false;

    private OwnColorTable ownColorTable; // Class to load LUT's

    /* Directionality stuff */


    @Override
    @SuppressWarnings("unchecked")
    public void run() {
        logService.info("Angle Analyzer 0.1");


        ownColorTable = new OwnColorTable(lutService);
        ownColorTable.setLut("physics.lut");


        if(imp == null)
            imp = WindowManager.getCurrentImage();

        //imp = ImageJFunctions.wrap((RandomAccessibleInterval<T>) dataset, "input");
        ImagePlus max_holder = imp.createImagePlus();

        if(imp.hasImageStack()){
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
            int width = imp.getWidth();
            int height = imp.getHeight();

            for(int y = 0; y < height; y += window){
                for(int x = 0; x < width; x += window){
                        ImagePlus imp_window = util.applyWindow(imp, x, y, width, height);
                        ImageProcessor allpeaks = util.obtainDistancesFFT(imp_window, 0.005f);
                        ImagePlus result = new ImagePlus("FFT", allpeaks);

                        int r_height = result.getHeight();
                        int r_width = result.getWidth();
                        result.setRoi(0, 0, r_width/2, r_height);
                        ProfilePlot profilePlot = new ProfilePlot(imp);
                        double[] profile = profilePlot.getProfile();

                        int index = 0;
                        for (int i = 0; i < profile.length; i++) {
                            index = profile[i] > profile[index] ? i : index;
                        }
                        float angle = index/((float)r_width/2)*180; //angle in degrees from 0-180

                        max_holder.setRoi(x, y, width, height);
                        Color color = ownColorTable.getColor(angle, 0f, 180f);
                        max_holder.setColor(color);
                        max_holder

                    //result.show();

                    }
                }
        }






        /*
        // threshold remove bottom 5%

        ImagePlus mask = imp.duplicate();
        util.MakeMask(mask);


        String options = "tensor=1.0 gradient=0 radian=on vectorgrid=10 vectorscale=300.0 vectortype=3 vectoroverlay=false vectortable=on";
        logService.info("Running OrientationJ using: " + options);
        if(debugging){
            logService.info("Running from IDE");
            WindowManager.setTempCurrentImage(imp);
            OrientationJ_Vector_Field vector_field = new OrientationJ_Vector_Field();

            //Macro.setOptions(options);
            vector_field.run(options);
        } else {
            IJ.run("OrientationJ Vector Field", options);
        }

        //This is bad and i should feel bad
        ResultsTable resultsTable = ResultsTable.getResultsTable("OJ-Table-Vector-Field-");

        System.out.println(resultsTable.getColumnHeadings());



         */

        System.out.println("done");
	} // substack 32-40

    // Only run from the IDE
    public static void main(final String... arguments) throws IOException {
        debugging = true;
        ImageJ ij = new ImageJ();
        ij.ui().showUI();

        //Dataset input = (Dataset) ij.io().open("test_stack.tif");
        //ij.ui().show(input);


        ImagePlus imp = new ij.io.Opener().openImage("W:\\Data\\Microscopy\\RCM\\Test Data\\SPC Horizontal\\MAX_middle to bottom- bottom 30% 405 5% 561_1_MMStack_Pos0.ome-1-1.tif");
        imp.show();
        //args = "D:\\Data\\Microscopy\\2022\\07\\8%561_40ms_MP3_1_RhB100x\\height slice_1\\slice_36_crop.tif";
        // "C:\\Users\\gobes001\\LocalSoftware\\AnalyseDirectionality\\test.tif";
        ij.command().run(AngleAnalyzer.class, true);
    }
}


