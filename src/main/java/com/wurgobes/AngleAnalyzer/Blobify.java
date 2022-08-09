package com.wurgobes.AngleAnalyzer;

import ij.ImagePlus;
import ij.ImageStack;
import ij.Macro;
import ij.WindowManager;
import ij.process.ImageProcessor;
import net.imagej.ImageJ;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

@Plugin(type = Command.class, name = "Blobify", menuPath = "Plugins>Blobify")
public class Blobify  implements Command {

    @Parameter
    private LogService logService;

    //@Parameter
    protected ImagePlus imp;

    private static String args;
    private static boolean debugging = false;
    private boolean headless = false;

    Blobify(ImagePlus imp) {
        this.imp = imp;
    }

    private void setup_image() {
        // null if no image is found
        imp = WindowManager.getCurrentImage();

        if(args.equals(""))
            args = Macro.getOptions();
        if(args != null && !args.equals(""))
            headless = true;

        if(imp == null) {
            imp = new ij.io.Opener().openImage(args);
            imp.show();
        }

        if(imp == null){
            logService.error("Expected Image");
            throw new IllegalArgumentException();
        }
        if(!isBinary(imp)) {
            logService.error("Expected Binary Image");
            throw new IllegalArgumentException();
        }
    }

    private boolean isBinary(ImagePlus imp) {
        ImageProcessor slice;
        if(imp.isStack()) {
            slice = imp.getStack().getProcessor(1);
        } else {
            slice = imp.getProcessor();
        }
        return slice.isBinary();
    }


    private void ProcessSlice(ImageProcessor imageProcessor){

    }
    @Override
    public void run() {
        setup_image();
        if(imp.isStack()){
            ImageStack imageStack = imp.getStack();
            int slices = imageStack.getSize();
            for(int n = 1; n <= slices; n++){
                ProcessSlice(imageStack.getProcessor(n));
            }

        }

    }

    // Only run from the IDE
    public static void main(final String... arguments) {
        debugging = true;
        final ImageJ ij = new ImageJ();

        ij.launch(arguments);

        args = "D:\\Data\\Microscopy\\Unilever\\crop_blob.tif";
        ij.command().run(Blobify.class, true);
    }
}
