package com.wurgobes.AngleAnalyzer;

import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.Thresholder;
import ij.process.ImageProcessor;
import net.imglib2.histogram.Histogram1d;
import net.imglib2.histogram.Integer1dBinMapper;

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
    public static void Mask(ImagePlus imp){
        if(imp.isStack()){
            ImageStack imageStack = imp.getStack();
            int slices = imageStack.getSize();
            for(int n = 1; n <= slices; n++){
                ImageProcessor imageProcessor = imageStack.getProcessor(n).duplicate();
                processSlice(imageProcessor);
            }
        } else {
            processSlice(imp.getProcessor().duplicate());
        }
    }

    private static void processSlice(ImageProcessor imageProcessor) {
        ImagePlus temp = new ImagePlus("temp", imageProcessor);
        temp.show();

        // Create a histogram out of the cropped image
        Histogram1d<T> hist = new Histogram1d<>(distanceCenter, new Integer1dBinMapper<>(0, 256, false));

        // Give us a min and max to mask the image with
        int minThreshold = getThresholdBin(0.99_86f, hist.toLongArray());
        int maxThreshold = 256;

        imageProcessor.setThreshold(threshold, Math.pow(2, imageProcessor.getBitDepth()), ImageProcessor.BLACK_AND_WHITE_LUT);
        (new Thresholder()).run("mask");
        // threshold



        // binary - close
        // iter 10 count 2
    }

}
