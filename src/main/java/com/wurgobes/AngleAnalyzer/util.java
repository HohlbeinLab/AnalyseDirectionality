package com.wurgobes.AngleAnalyzer;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.*;
import ij.plugin.filter.RankFilters;
import ij.process.AutoThresholder;
import ij.process.FHT;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import net.imglib2.Point;
import net.imglib2.util.Pair;
import net.imglib2.util.ValuePair;

import java.util.*;



public class util {

    private static int[] toInt(Object input, int bitdepth){
        switch(bitdepth){
            case 8:
                return toInt((byte[]) input);
            case 16:
                return (int[]) input;
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


        ImagePlus mask = new ImagePlus("temp", imageProcessor);



        imageProcessor.setAutoThreshold(AutoThresholder.Method.Huang, true, ImageProcessor.BLACK_AND_WHITE_LUT);

        mask.show();


    }

    public static ImageProcessor obtainDistancesFFT(ImagePlus imp, float angleprecision) {
        List<Pair<Double, List<Pair<Double, Double>>>> results =  new ArrayList<>();


        final FHT fht = new FHT();

        double radius = Math.min(imp.getHeight(), imp.getWidth())/2.;
        Point center = new Point(imp.getWidth()/2, imp.getHeight()/2);
        int sidelength = (int) (radius * Math.sqrt(2));
        double dx = imp.getCalibration().pixelHeight;
        final double[] freq = getX(lowerpower(sidelength), dx);
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
            //Plot plot = new Plot("FFT", "x", "y");
            //plot.add("line", freq, magnitudes);
            //plot.show();

            List<Pair<Double, Double>> peaks = filterpeaks(toDouble(magnitudes), freq, 0.5, 0.25, 1, 0.7, dx);
            if (peaks.size() > 0)
                results.add(new ValuePair<>(angle+(Math.PI/2 - 1), peaks));
        }
        imp.setRoi(getRoi(center, sidelength, end)); // Clear ugly line

        for(Pair<Double, List<Pair<Double, Double>>> angle: results){
            if(angle.getB().size() > 1) {
                IJ.log(String.format("Angle: %.2f deg (%.2f rad)", Math.toDegrees(angle.getA()), angle.getA()));

                List<Pair<Double, Double>> peaks = angle.getB();
                peaks.sort(Comparator.comparingDouble(Pair::getB));
                Collections.reverse(peaks);


                for (Pair<Double, Double> peak : peaks.subList(0, 2))
                    IJ.log(String.format("Found peak at %.2f nm with confidence %.3f (%.2f Hz)", 1000 / peak.getA(), peak.getB(), peak.getA()));
            }
        }


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
        int x0 = (int) (center.getIntPosition(0) + Math.cos((angle) * Math.PI) * sidelength/2);
        int y0 = (int) (center.getIntPosition(1) + Math.sin((angle) * Math.PI) * sidelength/2);
        int x1 = (int) (center.getIntPosition(0) + Math.cos((1+angle) * Math.PI) * sidelength/2);
        int y1 = (int) (center.getIntPosition(1) + Math.sin((1+angle) * Math.PI) * sidelength/2);

        return new RotatedRectRoi(x0, y0, x1, y1, sidelength);
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
}
