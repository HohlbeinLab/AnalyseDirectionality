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

import net.imagej.Dataset;
import net.imagej.DatasetService;
import net.imagej.ImageJ;

import fiji.analyze.directionality.Directionality_;

import net.imglib2.type.numeric.RealType;
import org.scijava.ItemIO;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.*;


@Plugin(type = Command.class, name = "Angle Analyzer", menuPath = "Plugins>Angle Analyzer<Analyze Angles")
public class AngleAnalyzer <T extends RealType<T>> implements Command {

    // The services are passed through from ImageJ automatically
    @Parameter
    private DatasetService datasetService;

    @Parameter
    private LogService logService;

    @Parameter(type = ItemIO.OUTPUT)
    private Dataset dataset;

    @Parameter(min = "0")
    private float angleStart = 0;

    @Parameter(max = "180")
    private float angleEnd = 180;

    @Parameter(min = "1")
    private int sobelSize = 5;

    // Debug parameters
    private static boolean runningFromIDE = false;
 

    @Override
    public void run() {
        System.out.println("test");
	}

    // Only run from the IDE
    public static void main(final String... args) throws Exception {
        runningFromIDE = true; //this is really dumb


        final ImageJ ij = new ImageJ();
        ij.launch(args);

        ij.command().run(AngleAnalyzer.class, true);
    }
}

