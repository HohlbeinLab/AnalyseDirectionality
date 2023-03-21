package com.wurgobes.AngleAnalyzer;


import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.*;
import ij.plugin.filter.RankFilters;
import ij.process.*;
import net.imglib2.Point;
import net.imglib2.util.Pair;
import net.imglib2.util.ValuePair;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.nio.file.Path;
import java.util.*;
import java.util.List;
import static ij.measure.Measurements.*;



public class util {

    private static int[] toInt(Object input, int bitdepth){
        switch(bitdepth){
            case 8:
                return toInt((byte[]) input);
            case 32:
                return toInt((long[]) input);
            default:
                return (int[]) input;
        }
    }

    private static int[] toInt(byte[] input) {
        int[] output = new int[input.length];
        for(int i = 0; i < input.length; i++)
            output[i] = input[i];
        return output;
    }
    private static int[] toInt(long[] input) {
        int[] output = new int[input.length];
        for(int i = 0; i < input.length; i++)
            output[i] = (int) input[i];
        return output;
    }

    public static long sum(int[] rx) {
        long sum = 0L;

        for (int l : rx)
            sum += l;

        return sum;
    }

    private static int getThresholdBin(float threshold, final int[] hist){
        // Get the amount of items that are under the threshold in the histogram
        int bins = hist.length;
        final long elementSum = sum(hist);

        threshold *= elementSum;

        long temp = 0L;
        for(int i = 0; i < bins; i++){
            temp += hist[i];
            if(temp > threshold) return i;
        }

        return 0;
    }

    /**
     * Takes an imageplus and masks away certain areas according to the imput parameters
     * NOTE: MODIFIES imp IN PLACE
     * */
    public static void MakeMask(ImagePlus imp){
        if(imp.isStack()){
            ImageStack imageStack = imp.getStack();
            int slices = imageStack.getSize();
            for(int n = 1; n <= slices; n++){
                processSlice(imageStack.getProcessor(n));
            }
        } else {
            processSlice(imp.getProcessor());
        }
    }

    private static void processSlice(ImageProcessor imageProcessor) {
        // min of 10px, max off 15px
        // min = 1, max = 2
        RankFilters rankFilters = new RankFilters();
        rankFilters.rank(imageProcessor, 10, 1);
        rankFilters.rank(imageProcessor, 15, 2);


        //ImagePlus mask = new ImagePlus("temp", imageProcessor);



        imageProcessor.setAutoThreshold(AutoThresholder.Method.Huang, true, ImageProcessor.BLACK_AND_WHITE_LUT);

        //mask.show();


    }

    public static ImageProcessor obtainDistancesFFT(ImagePlus imp, float angleprecision) {

        final FHT fht = new FHT();

        double radius = Math.min(imp.getHeight(), imp.getWidth())/2.;
        Point center = new Point(imp.getWidth()/2, imp.getHeight()/2);
        int sidelength = (int) (radius * Math.sqrt(2));
        double dx = imp.getCalibration().pixelHeight;
        final double[] freq = getX(lowerpower(sidelength), dx);

        // in pi radians
        int start = 0;
        int end = 2;
        int length = Math.abs(end-start);

        float[][] temp = new float[(int) Math.ceil(length/angleprecision)+1][freq.length];

        int index = 0;

        for(double angle = start; angle <= end; angle += angleprecision){
            RotatedRectRoi roi = getRoi(center, sidelength, angle);

            imp.setRoi(roi);

            float[] data = getProfile(imp, roi);

            float[] magnitudes = fht.fourier1D(data, FHT.HAMMING);
            temp[index++] = magnitudes;


        }
        imp.setRoi(getRoi(center, sidelength, end));

        //return results;
        return new FloatProcessor(temp);
    }

    /**
    Temporary solution because the RotatedRect Profileplot support is only in the daily built atm
     **/
    private static float[] getProfile(ImagePlus imp, RotatedRectRoi roi){
        double[] p = roi.getParams();
        Roi line = new Line(p[0], p[1], p[2], p[3]);
        line.setStrokeWidth(p[4]);
        line.setImage(imp);
        imp.setRoi(line, false);

        ProfilePlot profilePlot = new ProfilePlot(imp);
        return toFloat(profilePlot.getProfile());
    }

    public static float[] toFloat(double[] v){
        float[] result = new float[v.length];
        for(int i = 0; i < v.length; i++) result[i] = (float) v[i];
        return result;
    }

    public static double[] toDouble(float[] v){
        double[] result = new double[v.length];
        Arrays.setAll(result, i -> v[i]);
        return result;
    }

    public static double[] getX(int length, double dx){
        double[] x = new double[length];
        Arrays.setAll(x, i -> i/(2*length*dx));
        return x;
    }

    public static List<Pair<Double, Double>> filterpeaks(double[] magnitude, double[] freq, double min_freq, double freq_window, double conf_threshold, double conf_cutoff,double dx){
        List<Pair<Double, Double>> peaks = new ArrayList<>();
        if(freq_window > min_freq) min_freq = freq_window;
        final int di = (int) (freq_window*2*freq.length*dx);
        for(int i = (int) (min_freq*2*freq.length*dx); i < magnitude.length; i++){

            double window_mean = Arrays.stream(Arrays.copyOfRange(magnitude, i-di, i+di)).average().orElse(Double.NaN);
            double confidence = magnitude[i] / (conf_threshold * window_mean);
            if(confidence > conf_cutoff)
                peaks.add(new ValuePair<>(freq[i], confidence));
        }
        return peaks;
    }

    public static RotatedRectRoi getRoi(Point center, double sidelength, double angle){
        int x0 = (int) (center.getIntPosition(0) + Math.cos((1+angle) * Math.PI) * sidelength/2);
        int y0 = (int) (center.getIntPosition(1) + Math.sin((1+angle)* Math.PI) * sidelength/2);
        int x1 = (int) (center.getIntPosition(0) + Math.cos((0+angle) * Math.PI) * sidelength/2);
        int y1 = (int) (center.getIntPosition(1) + Math.sin((0+angle) * Math.PI) * sidelength/2);

        return new RotatedRectRoi(x0, y0, x1, y1, sidelength);
    }

    public static RotatedRectRoi getRoi(Point center, double sidelength, double angle, double sidelength2){
        int x0 = (int) (center.getIntPosition(0) + Math.cos((1+angle) * Math.PI) * sidelength/2);
        int y0 = (int) (center.getIntPosition(1) + Math.sin((1+angle)* Math.PI) * sidelength/2);
        int x1 = (int) (center.getIntPosition(0) + Math.cos((0+angle) * Math.PI) * sidelength/2);
        int y1 = (int) (center.getIntPosition(1) + Math.sin((0+angle) * Math.PI) * sidelength/2);

        return new RotatedRectRoi(x0, y0, x1, y1, sidelength2);
    }

    private static int lowerpower(int n){
        int power = 1;
        while(n > power){
            power *= 2;
        }
        return power/2;
    }

    public static ImagePlus applyWindow(ImagePlus imp, int x, int y, int size_x, int size_y){
        ImagePlus copy = imp.duplicate();
        return copy.crop(new Roi[]{new Roi(x, y, size_x, size_y)})[0];
    }

    public static void addLutLegend(ImageProcessor plot, OwnColorTable ct, String label, int width, double start, double end){
        int img_height = plot.getHeight();
        int img_width = plot.getWidth();
        plot.setLineWidth(2);
        // Add a LUT legend to the provided plot with the provided start, end and title
        for(int i = 0; i < width; i++){
            plot.setColor(ct.getColor(i, 0, width));
            plot.drawLine((int) (img_width*(0.01 + 0.0005 * (i+2))), (int) (img_height*0.93), (int) (img_width*(0.01 + 0.0005 * (i+2))), (int) (img_height*0.99));
        }
        plot.setColor(new Color(0, 0,0));

        int[] args = getAbsoluteCoords(new double[]{0.01,0.925, 0.01, 0.995}, img_width, img_height);
        plot.drawLine(args[0], args[1], args[2], args[3]);
        args = getAbsoluteCoords(new double[]{0.01 + 0.0005 * (width+2),0.925, 0.01 + 0.0005 * (width+2), 0.995}, img_width, img_height);
        plot.drawLine(args[0], args[1], args[2], args[3]);
        args = getAbsoluteCoords(new double[]{0.011, 0.99, 0.01 + 0.0005 * (width+2), 0.99}, img_width, img_height);
        plot.drawLine(args[0], args[1], args[2], args[3]);


        double diff = end - start;

        for(double i = start + 25; i < end; i += 25){
            double w = ((i-start)/diff)*width;
            args = getAbsoluteCoords(new double[]{0.01 + 0.0005 * (w+2), 0.986, 0.01 + 0.0005 * (w+2), 0.99}, img_width, img_height);
            plot.drawLine(args[0], args[1], args[2], args[3]);

        }

        for(double i = start + 50; i < end; i += 50){
            double w = ((i-start)/diff)*width;
            args = getAbsoluteCoords(new double[]{0.01 + 0.0005 * (w+2), 0.984, 0.01 + 0.0005 * (w+2), 0.99}, img_width, img_height);
            plot.drawLine(args[0], args[1], args[2], args[3]);
        }

        for(double i = start + 100; i < end; i += 100){
            double w = ((i-start)/diff)*width;
            args = getAbsoluteCoords(new double[]{0.01 + 0.0005 * (w+2), 0.982, 0.01 + 0.0005 * (w+2), 0.99}, img_width, img_height);
            plot.drawLine(args[0], args[1], args[2], args[3]);
        }

        plot.setFontSize(32);
        plot.setColor(new Color(255, 255, 255));
        plot.drawString(String.valueOf((int) start), (int) (0.001*img_width), (int) (0.92*img_height));
        plot.drawString(label, (int) (img_width*(0.001 + 0.00025 * (width-label.length()*2))), (int) (img_height*0.92));
        plot.drawString(String.valueOf((int) end), (int) (img_width*(0.001 + 0.0005 * (width+2))), (int) (0.92*img_height));
    }

    private static int[] getAbsoluteCoords(double[] coords, int width, int height){
        return new int[]{(int) (width * coords[0]), (int)(height * coords[1]), (int)(width * coords[2]), (int) (height * coords[3])};
    }

    public static void SaveCSV(final ArrayList<ArrayList<Double>> data, List<String> Headers, Path CSV_FILE_NAME) {
        try {
            FileWriter fileWrt = new FileWriter(CSV_FILE_NAME.toString());
            BufferedWriter bufferWrt = new BufferedWriter(fileWrt);
            bufferWrt.write(String.join(",", Headers) + "\n");

            StringBuilder s = new StringBuilder(100);
            for (int r = 0; r < data.size(); r++) {
                s.setLength(0);
                for (int c = 0; c < data.get(0).size(); c++) {
                    s.append(data.get(r).get(c));
                    if (c < data.get(0).size() - 1) {
                        s.append(",");
                    }
                }
                if (r <= data.size() - 1) {
                    s.append("\n");
                    bufferWrt.write(s.toString());
                }
            }

            bufferWrt.close();
        } catch (Exception e) {
            System.out.println("Could not save CSV.");
            e.printStackTrace();
        }
        System.gc();
    }

    public static double[] calculateStandardDeviation(ArrayList<Double> array) {

        // get the sum of array
        double sum = 0.0;
        for (double i : array) {
            sum += i;
        }

        // get the mean of array
        int length = array.size();
        double mean = sum / length;

        // calculate the standard deviation
        double standardDeviation = 0.0;
        for (double num : array) {
            standardDeviation += Math.pow(num - mean, 2);
        }

        return new double[]{mean, Math.sqrt(standardDeviation / length)};
    }

    public static ImageProcessor multiply(ImageProcessor imp1, ColorProcessor imp2, double factor){
        ColorProcessor result = new ColorProcessor(imp1.getWidth(), imp1.getHeight());
        double max_val = 1<<imp1.getBitDepth();
        for(int i = 0; i < imp1.getPixelCount(); i++){
            double grey_value = Math.pow(imp1.get(i)/max_val, factor);
            int c = imp2.get(i);

            int r = (int) (grey_value*((c&0xff0000)>>16));
            int g = (int) (grey_value*((c&0xff00)>>8));
            int b = (int) (grey_value*(c&0xff));

            result.set(i, (r<<16) + (g<<8) + b);
        }

        return result;
    }

    public static double getMedian(ImageProcessor imp, int x, int y, int width, int height){
        imp.setRoi(x, y, width, height);
        ImageStatistics stats = ImageStatistics.getStatistics(imp, AREA+MEAN+STD_DEV+MODE+MIN_MAX+RECT+MEDIAN, null);
        return stats.median;
    }

    public static ArrayList<Double> getValuesWindow(ArrayList<Double> lst, int index, int window_size){
        
    }

}
