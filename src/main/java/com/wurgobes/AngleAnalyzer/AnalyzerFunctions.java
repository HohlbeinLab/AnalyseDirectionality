package com.wurgobes.AngleAnalyzer;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.*;
import ij.process.*;
import net.imglib2.Point;
import org.apache.commons.lang3.tuple.Pair;


import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

import static com.wurgobes.AngleAnalyzer.util.*;
import static com.wurgobes.AngleAnalyzer.util.addLutLegend;

public class AnalyzerFunctions {

    static Pair<Integer, ImageStack> run(ImagePlus imp, ImagePlus max_imp, ImagePlus mask, ArrayList<ArrayList<Double>> csv_data, ArrayList<Double> angle_map, OwnColorTable circularLut, RFTParameters params){
        double min = Double.MAX_VALUE;
        double max = -1;

        int winspace = (int) Math.ceil(params.window*(1-params.overlap));
        int width_mod = params.window + params.buffer;
        int height_mod = params.window + params.buffer;
        int map_width = -1;
        ImageStack fft_stack = null;

        ColorProcessor max_ip = null;
        if(!params.scanning_range & !params.macro_mode)
            max_ip = (ColorProcessor) max_imp.getProcessor();

        for(int y = params.buffer; y <= params.height - height_mod; y += winspace){
            for(int x = params.buffer; x <= params.width - width_mod; x += winspace){

                ImagePlus imp_window = util.applyWindow(imp, x, y, params.window, params.window);
                ImageProcessor allpeaks = util.obtainDistancesFFT(imp_window, 0.005f);
                ImagePlus result = new ImagePlus("FFT", allpeaks);

                if(fft_stack==null)
                    fft_stack = new ImageStack(allpeaks.getWidth(), allpeaks.getHeight());

                fft_stack.addSlice(allpeaks);

                int r_height = result.getHeight();
                int r_width = result.getWidth();

                result.setRoi(0, 0, r_width, r_height);
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
                double angle = (double) index /(r_width)*180; //angle in degrees from 0-180

                if(!params.scanning_range & !params.macro_mode) {
                    Color color = circularLut.getColor(angle, 180f, 0);
                    max_ip.setColor(color);
                    max_ip.fillRect(x, y, params.window, params.window);
                    max_imp.repaintWindow();
                }


                double median = getMedian(mask.getProcessor(), x, y, params.window, params.window);

                angle_map.add(angle);
                ArrayList<Double> curr_data = new ArrayList<>(Arrays.asList((double) x, (double) y, (double) params.window, (double) params.window, (double) index, median, angle));
                curr_data.addAll(DoubleStream.of(profile).boxed().collect(Collectors.toCollection(ArrayList::new))); // profile data

                csv_data.add(curr_data);
            }
            if(map_width == -1) //size of one row of data
                map_width = csv_data.size();
        }
        return  Pair.of(map_width, fft_stack);
    }

    static ImageProcessor calcSigMap(ArrayList<ArrayList<Double>> csv_data, ArrayList<Boolean> sig_map, RFTParameters params) {
        ImageProcessor sig_ip = new ByteProcessor(params.width, params.height);
        //0, 1, 2,     3,      4,     5,           6,     7,         8+
        //x, y, width, height, index, mask median, angle, relevance, FT data
        for (ArrayList<Double> c : csv_data) {
            if (c.get(5) > (params.intensity_cutoff * 255)) {
                sig_map.add(Boolean.TRUE);
                sig_ip.setColor(Color.white);
            } else {
                sig_map.add(Boolean.FALSE);
                sig_ip.setColor(Color.black);
            }

            int x = c.get(0).intValue();
            int y = c.get(1).intValue();
            sig_ip.fillRect(x, y, params.window, params.window);

        }
        return sig_ip;
    }

    static void calcVectorMap(ArrayList<Roi> rois, ArrayList<ArrayList<Double>> csv_data, RFTParameters params){
        //0, 1, 2,     3,      4,     5,           6,     7,         8+
        //x, y, width, height, index, mask median, angle, relevance, FT data
        for (ArrayList<Double> c : csv_data) {
            double ad = c.get(7);

            if (c.get(7) > params.cutoff && c.get(5) > (params.intensity_cutoff * 255)) {  // Significant
                Line line = getLineRoi(new Point((long) (c.get(0) + c.get(2) / 2), (long) (c.get(1) + c.get(3) / 2)), 0.1 * ad * params.window * params.vector_length, c.get(6) / 180f);
                line.setStrokeWidth(params.vector_width);
                rois.add(line);
            }
        }
    }

    static void applyOverlay(ImagePlus imp, Overlay overlay){
        imp.setOverlay(overlay);
    }

    static void applyOverlay(ImagePlus imp, Overlay overlay, ArrayList<Roi> rois, RFTParameters params){
        overlay.clear();
        overlay.setStrokeColor(params.vector_color);
        for(Roi roi : rois)
            overlay.add(roi);
        imp.setOverlay(overlay);

    }

    static ImageStack order_parameter(ArrayList<ArrayList<Double>> csv_data, ArrayList<Double> angle_map, ArrayList<Boolean> sig_map, OwnColorTable orderLUT, Plot order_plot, int map_width, int width, int height, int window){
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

        ArrayList<Double> angle_map_filtered = maskList(angle_map, sig_map, Double.NaN);

        ArrayList<Double> average_order =  new ArrayList<>();
        ArrayList<Integer> neigh_sizes =  new ArrayList<>();


        for(int neighbourhood_size = 3; neighbourhood_size < Math.min(map_length, map_width); neighbourhood_size += 2) {
            ArrayList<Double> order_parameter_submap = new ArrayList<>(angle_map.size());
            neigh_sizes.add(neighbourhood_size);

            int half_window = neighbourhood_size % 2 == 1 ? (neighbourhood_size-1)/2 : neighbourhood_size/2;
            for(int i = 0; i < angle_map.size(); i++) {
                // skip edge values
                if(
                        i % map_width < half_window || // left
                        i % map_width >= map_width - half_window || // right
                        i < map_width * half_window || // top
                        i / map_width > map_length - half_window // bottom
                )
                    order_parameter_submap.add(Double.NaN);
                else {
                    ArrayList<Double> value_window = getValuesWindow(angle_map_filtered, i, neighbourhood_size, map_width);
                    double central_value = value_window.remove(value_window.size()/2); // Remove Central Value
                    if(Double.isNaN(central_value)){
                        order_parameter_submap.add(Double.NaN);
                    } else {
                        ArrayList<Double> cos2window = ListCos(value_window, central_value);

                        double order_param = 2*(mean(cos2window)-0.5);
                        order_parameter_submap.add(order_param);
                    }
                }
            }
            average_order.add(mean(order_parameter_submap));
            order_parameter_map.add(order_parameter_submap);
        }


        ImageStack order_stack = new ImageStack(width, height);



        int neighbourhood_size = 3;
        for (int i = 0; i < order_parameter_map.size(); i++) { // single neighbourhood size
            // "x", "y", "width", "height", "Max Index", "Mask Median", "Angle", "Relevance?", "Profile Data"
            ArrayList<Double> order_submap = order_parameter_map.get(i);
            ColorProcessor order_ip = new ColorProcessor(width, height);

            for (int j = 0; j < order_submap.size(); j++) { // individual x_y
                ArrayList<Double> metadata = csv_data.get(j);
                int x = metadata.get(0).intValue();
                int y = metadata.get(1).intValue();

                double value = order_submap.get(j);
                Color color;
                if (Double.isNaN(value)) {
                    color = Color.black;
                } else {
                    color = orderLUT.getColor(order_submap.get(j), 0.0, 1.0);
                }

                order_ip.setColor(color);
                order_ip.fillRect(x, y, window, window);
            }
            addLutLegend(order_ip, orderLUT, "Order Parameter", 1024, 0.0, 1.0);

            order_stack.addSlice("NH " + neighbourhood_size, order_ip, i);
            neighbourhood_size += 2;

        }
        float[] average_order_hist = new float[average_order.size()];
        for(int i = 0; i < average_order.size(); i++)
            average_order_hist[i] =  average_order.get(i).floatValue();

        float[] average_order_hist_x = new float[neigh_sizes.size()];
        for(int i = 0; i < neigh_sizes.size(); i++)
            average_order_hist_x[i] =  neigh_sizes.get(i).floatValue();


        order_plot.addPoints(average_order_hist_x, average_order_hist, null, Plot.toShape("line"), "");


        return order_stack;
    }

    static void calc_adjusted_stats(ArrayList<Double> adjusted_stats, ArrayList<ArrayList<Double>> csv_data){
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

        //0, 1, 2,     3,      4,     5,           6,     7,         8+
        //x, y, width, height, index, mask median, angle, relevance, FT data
        for (int i = 0; i < csv_data.size(); i++) {
            ArrayList<Double> c = csv_data.get(i);
            double ad = adjusted_stats.get(i);
            c.add(7, ad);
            csv_data.set(i, c);
        }
    }

    static void createAngleGraph(ImageStack fft_stack, ArrayList<ArrayList<Double>> csv_data, double cutoff, Plot hist_alt){
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
            fft_dummy.setRoi(0, 3, (fft_dummy.getWidth()/2), Math.min(fft_dummy.getHeight(), 22));
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
        hist_alt.setColor(Color.red);
        hist_alt.addPoints(xValues_fft, histogram_fft_o, null, Plot.toShape("line"), "Over Threshold");
        hist_alt.setColor(Color.blue);
        hist_alt.addPoints(xValues_fft, histogram_fft_u, null, Plot.toShape("line"), "Under Threshold");
        hist_alt.setColor(Color.green);
        hist_alt.addPoints(xValues_fft, histogram_fft_a, null, Plot.toShape("line"), "Total");

        hist_alt.addLegend(null);
        hist_alt.setLimitsToFit(true);
    }

    static void createMultImage(ImagePlus imp, ImagePlus mask, ColorProcessor max_ip, OwnColorTable circularLut){
        ImagePlus imp_clone = imp.duplicate();
        ImageConverter converter = new ImageConverter(imp_clone);
        converter.convertToRGB();
        ImageProcessor mul_ip = multiply(mask.getProcessor(), max_ip, 2);
        addLutLegend(mul_ip, circularLut, "Angle", 1024, 0f, 180f);
        ImagePlus mul_result = new ImagePlus("mul_result", mul_ip);
        mul_result.show();
    }


    static void toggle_overlay(ImagePlus imp, Overlay overlay){
        //Doesnt work
        if(imp.getOverlay().size() == 0) {
            applyOverlay(imp, overlay);
        } else {
            imp.setOverlay(new Overlay());
        }
    }
}
