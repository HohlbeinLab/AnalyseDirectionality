package com.wurgobes.RFT;

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

import static com.wurgobes.RFT.util.*;
import static com.wurgobes.RFT.util.addLutLegend;

public class AnalyzerFunctions {

    static Pair<Integer, ImageStack> run(ImagePlus imp, ImagePlus max_imp, ImagePlus mask, ArrayList<ArrayList<Double>> csv_data, ArrayList<Double> angle_map, OwnColorTable circularLut, RFTParameters params){


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
                double[] profile = new ProfilePlot(result).getProfile();

                int index = 0;
                for (int i = 0; i < profile.length; i++) {
                    index = profile[i] > profile[index] ? i : index;
                }
                double peak_angle = (double) index /(r_width)*180; //angle in degrees from 0-180

                if(!params.scanning_range & !params.macro_mode) {
                    Color color = circularLut.getColor(peak_angle, 180f, 0);
                    max_ip.setColor(color);
                    max_ip.fillRect(x, y, params.window, params.window);
                    max_imp.repaintWindow();
                }


                double median = getMedian(mask.getProcessor(), x, y, params.window, params.window);

                angle_map.add(peak_angle);
                ArrayList<Double> curr_data = new ArrayList<>(Arrays.asList((double) x, (double) y, (double) params.window, (double) params.window, (double) index, median, peak_angle));
                curr_data.addAll(DoubleStream.of(profile).boxed().collect(Collectors.toCollection(ArrayList::new))); // profile data

                csv_data.add(curr_data);
            }
            if(map_width == -1) //size of one row of data
                map_width = csv_data.size();
        }
        return  Pair.of(map_width, fft_stack);
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
        int hist_size = fft_stack.getWidth();
        float[] histogram_fft_o = new float[hist_size];
        float[] histogram_fft_u = new float[hist_size];
        float[] histogram_fft_a = new float[hist_size];
        float[] xValues_fft = new float[hist_size];
        for(int i = 0; i < xValues_fft.length; i++)
            xValues_fft[i]=(i/(float)xValues_fft.length)*180;

        for(int i = 1; i <= fft_stack.size(); i++){
            ImageProcessor fft_slice = fft_stack.getProcessor(i);
            ImagePlus fft_dummy = new ImagePlus("dummy", fft_slice);
            fft_dummy.setRoi(0, 3, fft_dummy.getWidth(), Math.min(fft_dummy.getHeight(), 22));
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

        hist_alt.setColor(Color.green);
        hist_alt.addPoints(xValues_fft, histogram_fft_a, null, Plot.toShape("line"), "Total");
        hist_alt.setColor(Color.red);
        hist_alt.addPoints(xValues_fft, histogram_fft_o, null, Plot.toShape("line"), "Over Threshold");
        hist_alt.setColor(Color.blue);
        hist_alt.addPoints(xValues_fft, histogram_fft_u, null, Plot.toShape("line"), "Under Threshold");


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
