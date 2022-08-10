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
import ij.Macro;
import ij.WindowManager;
import ij.measure.ResultsTable;
import net.imagej.ImageJ;
import net.imglib2.type.numeric.RealType;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import OrientationJ.OrientationJ_Vector_Field;


@Plugin(type = Command.class, name = "Angle Analyzer", menuPath = "Plugins>Angle Analyzer>Analyze Angles")
public class AngleAnalyzer <T extends RealType<T>> implements Command {

    // The services are passed through from ImageJ automatically

    @Parameter
    private LogService logService;



    /** The ImagePlus this plugin operates on. */
    //@Parameter(type= ItemIO.INPUT)
    protected ImagePlus imp;


    /* private variables */

    private static String args = "";

    private boolean headless = false;

    private static boolean debugging = false;

    /* Directionality stuff */


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
    }


    @Override
    public void run() {
        logService.info("Angle Analyzer 0.1");

        setup_image();

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




	} // substack 32-40

    // Only run from the IDE
    public static void main(final String... arguments) {
        debugging = true;
        final ImageJ ij = new ImageJ();

        ij.launch(arguments);

        //args = "D:\\Data\\Microscopy\\2022\\07\\8%561_40ms_MP3_1_RhB100x\\height slice_1\\slice_36_crop.tif";
        args = "E:\\Data\\Microscopy\\RCM\\2022\\07\\8%561_40ms_MP3_1_RhB100x\\height slice_1\\slice_xx.tif";
        ij.command().run(AngleAnalyzer.class, true);
    }
}


