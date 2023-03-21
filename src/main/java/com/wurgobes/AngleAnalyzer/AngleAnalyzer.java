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
import ij.measure.Measurements;
import ij.process.ColorProcessor;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import net.imagej.ImageJ;
import net.imagej.lut.LUTService;
import net.imglib2.Point;
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

    @Parameter(label = "Overlap", min = "0", max = "1", stepSize = "0.05")
    private double overlap = 0.5; // percent overlap

    @Parameter(label = "Buffer")
    private double buffer = 1; // percent overlap

    @Parameter(label = "Length", min = "0", max = "10", stepSize = "0.1")
    private double vector_length = 0.5;

    @Parameter(label = "Thickness", min = "0", max = "100", stepSize = "0.05")
    private double vector_thickness = 0.05;

    @Parameter(label = "Over/Under Cutoff (std over mean)", min = "0", max="5", stepSize ="0.1")
    private double cutoff = 3;

    @Parameter(label = "Intensity cutoff (0-255)", min = "0", max = "255", stepSize ="1")
    private double intensity_cutoff = 0;


    private Color vector_color = new Color(255, 227, 0);


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
        ColorProcessor max_ip = new ColorProcessor(width, height);
        ImagePlus max_imp = new ImagePlus("Colored blocks", max_ip);
        max_imp.show();

        // Mask
        ImagePlus mask = imp.duplicate();
        mask.setTitle("Mask");
        util.MakeMask(mask);
        ImageConverter converter = new ImageConverter(mask);
        converter.convertToGray8();
        mask.show();

        ImageStack fft_stack = null;

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


            int width_mod = (int) (width%(window*overlap));
            int height_mod = (int) (height%(window*overlap));

            double iter_max = window*Math.max(1-overlap + buffer, 1);

            double min = Double.MAX_VALUE;
            double max = -1;

            Overlay overlay = new Overlay();
            imp.setOverlay(overlay);
            //x, y, width, height, index, median, angle, FT data
            ArrayList<ArrayList<Double>> csv_data = new ArrayList<>();

            ArrayList<Double> angle_map = new ArrayList<>();
            int map_width = -1;

            logService.info("window: " + window +", width_mod: " + width_mod + ", height_mod: " + height_mod + ", buffer: " + buffer + ", iter_max: " + iter_max + ", overlap: " + overlap);
            for(int y = height_mod/2; y <= height - iter_max; y += (window*(1-overlap))){
                for(int x = width_mod/2; x <= width - iter_max; x += (window*(1-overlap))){


                        ImagePlus imp_window = util.applyWindow(imp, x, y, window, window);
                        ImageProcessor allpeaks = util.obtainDistancesFFT(imp_window, 0.005f);
                        ImagePlus result = new ImagePlus("FFT", allpeaks);

                        if(fft_stack==null)
                            fft_stack = new ImageStack(allpeaks.getWidth(), allpeaks.getHeight());

                        fft_stack.addSlice(allpeaks);

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


                    Color color = ownColorTable.getColor(angle, 180f, 0);
                    max_ip.setColor(color);
                    max_ip.fillRect(x, y, window, window);
                    max_imp.repaintWindow();

                    double median = getMedian(mask.getProcessor(), x, y, window, window);

                    angle_map.add((double) angle);
                    ArrayList<Double> curr_data = new ArrayList<>(Arrays.asList((double) x, (double) y, (double) window, (double) window, (double) index, median, (double) angle));
                    ArrayList<Double> profile_data = DoubleStream.of(profile).boxed().collect(Collectors.toCollection(ArrayList::new));
                    curr_data.addAll(profile_data); // java bad
                    curr_data.addAll(profile_data); // java bad

                    csv_data.add(curr_data);

                    //result.show();

                }
                if(map_width == -1) //size of one row of data
                    map_width = csv_data.size();
            }


            // calculate significance
            ArrayList<Double> adjusted_stats = new ArrayList<>();


            for (ArrayList<Double> list : csv_data){
                ArrayList<Double> adjusted_list = new ArrayList<>();
                for(int i = 7; i < list.size(); i++)
                    adjusted_list.add(list.get(i));

                double[] stats = calculateStandardDeviation(adjusted_list); //mean, std
                double cur_max = Collections.max(adjusted_list);

                if(Double.isNaN(stats[1]) || stats[1] == 0.0)
                    adjusted_stats.add(0.0);
                else {
                    double sig = ((cur_max - stats[0]) / stats[1]); //standard deviations over mean
                    adjusted_stats.add(sig);
                }
            }

            ArrayList<Boolean> sig_map = new ArrayList<>();


            //0, 1, 2,     3,      4,     5,           6,     7,         8+
            //x, y, width, height, index, mask median, angle, relevance, FT data
            for (int i = 0; i < csv_data.size(); i++){


                ArrayList<Double> c = csv_data.get(i);
                double ad = adjusted_stats.get(i);
                c.add(7, ad);
                csv_data.set(i, c);

                if(c.get(7) > cutoff)  // Significant
                    overlay.add(getRoi(new Point((long) (c.get(0) + c.get(2) / 2), (long) (c.get(1) + c.get(3) / 2)), 0.25 * window * vector_thickness, c.get(6) / 180f, 0.1 * ad * window * vector_length));
                if(c.get(5) > intensity_cutoff) {
                    sig_map.add(Boolean.TRUE);
                } else {
                    sig_map.add(Boolean.FALSE);
                }

            }



            //FFT stuff
            assert fft_stack != null;
            int hist_size = fft_stack.getWidth()/2;
            float[] histogram_fft_o = new float[hist_size];
            float[] histogram_fft_u = new float[hist_size];
            float[] histogram_fft_a = new float[hist_size];
            float[] xValues_fft = new float[hist_size];
            for(int i = 0; i < xValues_fft.length; i++)
                xValues_fft[i]=(i/(float)xValues_fft.length)*180;

            for(int i = 1; i <= fft_stack.size(); i++){
                ImageProcessor fft_slice = fft_stack.getProcessor(i);
                ImagePlus fft_dummy = new ImagePlus("dummy", fft_slice);
                fft_dummy.setRoi(0, 3, (fft_dummy.getWidth()/2), 14);
                ProfilePlot profilePlot = new ProfilePlot(fft_dummy);
                double[] profile = profilePlot.getProfile();

                if (csv_data.get(i-1).get(7) > cutoff)
                    for(int j = 0; j < profile.length; j++) {
                        histogram_fft_o[profile.length - j - 1] += profile[j];
                        histogram_fft_a[profile.length - j - 1] += profile[j];
                    }
                else
                    for(int j = 0; j < profile.length; j++) {
                        histogram_fft_u[profile.length - j - 1] += profile[j];
                        histogram_fft_a[profile.length - j - 1] += profile[j];
                    }
            }



            Plot hist_alt = new Plot("Histogram", "Angle", "Value");
            hist_alt.setColor(Color.red);
            hist_alt.addPoints(xValues_fft, histogram_fft_o, null, Plot.toShape("line"), "Over Threshold");
            hist_alt.setColor(Color.blue);
            hist_alt.addPoints(xValues_fft, histogram_fft_u, null, Plot.toShape("line"), "Under Threshold");
            hist_alt.setColor(Color.green);
            hist_alt.addPoints(xValues_fft, histogram_fft_a, null, Plot.toShape("line"), "Total");

            hist_alt.addLegend(null);
            hist_alt.setLimitsToFit(true);
            hist_alt.show();

            ImagePlus fft_imp = new ImagePlus("ffts", fft_stack);
            fft_imp.show();

            IJ.saveAsTiff(fft_imp, ".\\fft_stack.tif");

            converter = new ImageConverter(imp);
            converter.convertToRGB();
            overlay.fill(imp, vector_color, null);
            overlay.clear();
            SaveCSV(csv_data, new ArrayList<>(Arrays.asList("x", "y", "width", "height", "Max Index", "Mask Median", "Angle", "Relevance?", "Profile Data")), Paths.get(".\\test.csv"));

            /*
             * neighbourhood size search
             * from min to max size, steps of 2
             * heavily depends on window size
             * max size is maximum square that can fit (min of length and width)
             *
             * map_width - width of map
             * map_length (angle_map.size()/map_width) - length of map
             */

            int map_length = angle_map.size()/map_width;
            ArrayList<ArrayList<Double>> order_parameter_map = new ArrayList<>(); // per neighbourhood size (3, 5, 7, ...) a flat list of order parameters
            for(int neighbourhood_size = 3; neighbourhood_size < Math.min(map_length, map_width); neighbourhood_size += 2) {
                ArrayList<Double> order_parameter_submap = new ArrayList<>();

                for(int i = 0; i < angle_map.size(); i++) {

                }
            }

            ImageProcessor mul_ip = multiply(mask.getProcessor(), max_ip, 2);

            addLutLegend(mul_ip, ownColorTable, "Angle", 1024, 0f, 180f);
            addLutLegend(max_ip, ownColorTable, "Angle", 1024, 0f, 180f);
            ImagePlus mul_result = new ImagePlus("mul_result", mul_ip);
            mul_result.show();
            max_imp.repaintWindow();

        }









        /*
        // threshold remove bottom 5%




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
       // ImagePlus imp = new ij.io.Opener().openImage("W:\\Data\\Processing Testing\\MRI\\anisotropic_cropped.tif");
        //ImagePlus imp = new ij.io.Opener().openImage("W:\\Data\\Processing Testing\\Test Images\\tif\\45degright.tif");
        ImagePlus imp = new ij.io.Opener().openImage("test.tif");
        //ImagePlus imp = new ij.io.Opener().openImage("W:\\Data\\Microscopy\\Airyscan\\2022\\12\\Dead Stop SPC June 2022 11 cm\\10_crop.tif");
        //ImagePlus imp = new ij.io.Opener().openImage("W:\\Data\\Microscopy\\Airyscan\\2022\\12\\Dead Stop SPC June 2022 24 cm\\23_crop.tif");
        //ImagePlus imp = new ij.io.Opener().openImage("W:\\Data\\Microscopy\\Airyscan\\TVP_1\\TVP.tif");
        imp.show();
        //args = "D:\\Data\\Microscopy\\2022\\07\\8%561_40ms_MP3_1_RhB100x\\height slice_1\\slice_36_crop.tif";
        // "C:\\Users\\gobes001\\LocalSoftware\\AnalyseDirectionality\\test.tif";
        ij.command().run(AngleAnalyzer.class, true);
    }
}


