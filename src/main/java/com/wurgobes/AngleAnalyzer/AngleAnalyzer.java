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
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.ProfilePlot;
import ij.process.ImageConverter;
import net.imglib2.Point;
import ij.process.ColorProcessor;
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
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

import static com.wurgobes.AngleAnalyzer.util.*;


@Plugin(type = Command.class, name = "Angle Analyzer", menuPath = "Plugins>Angle Analyzer>Analyze Angles", priority = Priority.HIGH)
public class AngleAnalyzer <T extends RealType<T>> implements Command {

    // The services are passed through from ImageJ automatically

    @Parameter
    private LogService logService;

    @Parameter
    private LUTService lutService;

    @Parameter(label = "Size of processing window")
    private int window = 300; //px

    @Parameter(label = "Overlap")
    private double overlap = 0.5; // percent overlap

    @Parameter(label = "Buffer")
    private double buffer = 1; // percent overlap

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

        int width = imp.getWidth();
        int height = imp.getHeight();
        //imp = ImageJFunctions.wrap((RandomAccessibleInterval<T>) dataset, "input");
        //Non stack version for now
        ImageProcessor max_ip = new ColorProcessor(width, height);
        ImagePlus max_holder = new ImagePlus("Dominant Direction",  max_ip);
        max_holder.show();

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

            int totals = (height/window) * (width/window);
            int width_mod = (int) Math.max(width % window, buffer*window);
            int height_mod = (int) Math.max(height % window, buffer*window);

            double iter_max = window*Math.max(1-overlap + buffer, 1);

            int cnt = 0;

            double min = Double.MAX_VALUE;
            double max = -1;

            Overlay overlay = new Overlay();
            imp.setOverlay(overlay);
            //x, y, width, height, index, angle, FT data
            ArrayList<ArrayList<Double>> csv_data = new ArrayList<>();

            for(int y = height_mod/2; y <= height- iter_max; y += (window*(1-overlap))){
                for(int x = width_mod/2; x <= width - iter_max; x += (window*(1-overlap))){
                        cnt++;
                        logService.info(cnt/(double)totals);
                        ImagePlus imp_window = util.applyWindow(imp, x, y, window, window);
                        ImageProcessor allpeaks = util.obtainDistancesFFT(imp_window, 0.005f);
                        ImagePlus result = new ImagePlus("FFT", allpeaks);
                        int r_height = result.getHeight();
                        int r_width = result.getWidth();
                        result.setRoi(0, 0, r_width/2, r_height);
                        ProfilePlot profilePlot = new ProfilePlot(result);
                        double[] profile = profilePlot.getProfile();

                        int index = 0;
                        for (int i = 0; i < profile.length; i++) {
                            index = profile[i] > profile[index] ? i : index;
                            if (profile[i] < min)
                                min = profile[i];
                            else if (profile[i] > max)
                                max = profile[i];

                        }
                    float angle = index/(r_width/2f)*180; //angle in degrees from 0-180

                    angle = (angle - 90)*-1 + 90;

                    Color color = ownColorTable.getColor(angle, 0f, 180f);
                    max_ip.setColor(color);
                    max_ip.fillRect(x, y, window, window);
                    max_holder.repaintWindow();

                    ArrayList<Double> curr_data = new ArrayList<>(Arrays.asList((double) x, (double) y, (double) window, (double) window, (double) index, (double) angle));
                    ArrayList<Double> profile_data = DoubleStream.of(profile).boxed().collect(Collectors.toCollection(ArrayList::new));
                    curr_data.addAll(profile_data); // java bad
                    curr_data.addAll(profile_data); // java bad

                    csv_data.add(curr_data);

                    //result.show();

                }
            }

            ArrayList<Double> adjusted_stats = new ArrayList<>();

            for (ArrayList<Double> list : csv_data){;
                ArrayList<Double> adjusted_list = new ArrayList<>();
                for(int i = 6; i < list.size(); i++)
                    adjusted_list.add((list.get(i)-min)/(max-min));

                double[] stats = calculateStandardDeviation(adjusted_list); //mean, std
                double cur_max = Collections.max(adjusted_list);
                if(Double.isNaN(stats[1]) || stats[1] == 0.0)
                    adjusted_stats.add(0.0);
                else
                    adjusted_stats.add(((cur_max-stats[0])/stats[1])*cur_max); //(max-mean)/std
            }

            double stat_max = Collections.max(adjusted_stats);


            int binwidth = 5;
            double[] xValues = new double[180/binwidth + 1];
            for(int i = 0; i < xValues.length; i++)
                xValues[i]=i*binwidth;

            double[] histogram = new double[181];

            for (int i = 0; i < csv_data.size(); i++){
                double ad = adjusted_stats.get(i)/stat_max;
                ArrayList<Double> c = csv_data.get(i);
                c.add(6, ad);
                csv_data.set(i, c);
                overlay.add(getRoi(new Point((long) (c.get(0)+c.get(2)/4), (long) (c.get(1)+c.get(3)/4)), ad*window/2f, c.get(5)/180f, window/25f));

                histogram[(int)Math.round(c.get(5))/binwidth] += ad;
            }
            Plot plot = new Plot("Histogram", "Angle", "Value");
            plot.add("line", xValues, histogram);
            plot.show();

            ImageConverter converter = new ImageConverter(imp);
            converter.convertToRGB();
            overlay.fill(imp, new Color(255, 227, 0), null);
            SaveCSV(csv_data, new ArrayList<>(Arrays.asList("x", "y", "width", "height", "Max Index", "Angle", "Relevance?", "Profile Data")), Paths.get("W:\\Data\\Processing Testing\\test.csv"));
        }
        addLutLegend(max_ip, ownColorTable, "Angle", 1024, 0f, 180f);
        max_holder.setProcessor(max_ip);
        max_holder.repaintWindow();


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


        //ImagePlus imp = new ij.io.Opener().openImage("W:\\Data\\Microscopy\\RCM\\Test Data\\SPC Horizontal\\MAX_middle to bottom- bottom 30% 405 5% 561_1_MMStack_Pos0.ome-1-1.tif");
        ImagePlus imp = new ij.io.Opener().openImage("W:\\Data\\Processing Testing\\Test Images\\tif\\mosaic2.tif");
        //ImagePlus imp = new ij.io.Opener().openImage("W:\\Data\\Microscopy\\Airyscan\\2022\\12\\Dead Stop SPC June 2022 11 cm\\17.tif");
        imp.show();
        //args = "D:\\Data\\Microscopy\\2022\\07\\8%561_40ms_MP3_1_RhB100x\\height slice_1\\slice_36_crop.tif";
        // "C:\\Users\\gobes001\\LocalSoftware\\AnalyseDirectionality\\test.tif";
        ij.command().run(AngleAnalyzer.class, true);
    }
}


