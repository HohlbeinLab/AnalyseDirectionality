package com.wurgobes.AngleAnalyzer;

import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.filter.RankFilters;
import ij.process.AutoThresholder;
import ij.process.ImageProcessor;


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

}
