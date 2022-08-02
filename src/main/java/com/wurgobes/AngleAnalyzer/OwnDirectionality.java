package com.wurgobes.AngleAnalyzer;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.Point;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.util.ArrayList;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;

import fiji.analyze.directionality.Directionality_;
import ij.*;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.LookupPaintScale;
import org.jfree.chart.renderer.xy.ClusteredXYBarRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import ij.gui.ImageCanvas;
import ij.gui.Line;
import ij.gui.NewImage;
import ij.gui.Roi;
import ij.gui.TextRoi;

import ij.measure.CurveFitter;
import ij.measure.ResultsTable;

import ij.plugin.Duplicator;
import ij.plugin.filter.Convolver;
import ij.plugin.filter.GaussianBlur;

import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FHT;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

@Plugin(type = Command.class, name = "OwnDirectionality", menuPath = "Analyze>OwnDirectionality")
public class OwnDirectionality implements Command
{
    @Parameter
    private LogService logService;

    /*
     * ENUMS
     */

    public enum AnalysisMethod
    {
        FOURIER_COMPONENTS,
        LOCAL_GRADIENT_ORIENTATION;
        @Override
        public String toString()
        {
            switch ( this )
            {
                case FOURIER_COMPONENTS:
                    return "Fourier Components";
                case LOCAL_GRADIENT_ORIENTATION:
                    return "Local Gradient Orientation";
            }
            return "Not implemented";
        }

        public String toCommandName()
        {
            switch ( this )
            {
                case FOURIER_COMPONENTS:
                    return "Fourier";
                case LOCAL_GRADIENT_ORIENTATION:
                    return "Gradient";
            }
            return "Not implemented";
        }
    }

    /*
     * FIELDS
     */

    /* CONSTANTS */

    // Get rid of pixels too close to the center in the FFT spectrum.
    private static final float FREQ_THRESHOLD = 5.0f;

    /**
     * How many sigmas away from the gaussian center we sum to get the amount
     * value.
     */
    private static final double SIGMA_NUMBER = 2;

    private static final String PLUGIN_NAME = "Directionality analysis";

    private static final String VERSION_STR = "3.0.0";

    /* USER SETTING FIELDS, memorized between runs */

    /** The ImagePlus this plugin operates on. */
    @Parameter(label = "Image", description = "The ImagePlus this plugin operates on")
    protected ImagePlus imp;

    /** Method used for analysis, as set by the user. */
    @Parameter(choices = {"Fourier Components", "Local Gradient Orientation"}, label = "Method", description = "Method used for analysis")
    private static String string_method = "Fourier Components";

    /** The number of bins to create. */
    @Parameter(min = "1", label = "NBins", description = "The number of bins to create", stepSize = "1")
    private static int nbins = 90;

    /**
     * The first bin in degrees when displaying the histogram, so that we are
     * not forced to start at -90
     */
    @Parameter(label = "Histogram Start", description = "The first bin in degrees when displaying the histogram", stepSize = "1")
    private static double bin_start = -90;

    /**
     * The last bin in degrees when displaying the histogram, so that we are
     * able to limit the range of angles to analyse
     */
    @Parameter(label = "Histogram End", description = "The last bin in degrees when displaying the histogram", stepSize = "1")
    private static double bin_end = 90;

    @Parameter(min = "3", label = "Sobel Size", description = "The size of the sobel filter to use", stepSize = "2")
    private static int sobel_size = 5;

    /** If true, will calculate a map of orientation. */
    @Parameter(label = "Display orientation image", description = "Calculate an image with colors showing angle")
    private static boolean build_orientation_map = false;

    /** If true, will display a color wheel to interpret the orientation map. */
    @Parameter(label = "Display color wheel", description = "Show a wheel showing which angle is which color")
    private static boolean display_color_wheel = false;

    /**  If set true, will display a {@link ResultsTable} with the histogram at the end of processing. */
    @Parameter(label = "Display histogram table", description = "Display a table showing histogram statistics")
    private static boolean display_table = false;

    /** Debugging settings */
    @Parameter(label = "Debug", description = "If true, will display FFTs and filters")
    private static boolean debug = false;

    /* STD FIELDS */

    /** Enum to keep track of type. */
    private AnalysisMethod method = AnalysisMethod.FOURIER_COMPONENTS;

    /** FloatProcessor to convert source ImageProcessor to. */
    private FloatProcessor fip;

    /** Fourier filters are stored as a stack */
    protected ImageStack filters;

    /** Polar coordinates, stored as a FloatProcessor. */
    protected FloatProcessor window, r, theta;

    protected int width, height, small_side, long_side, npady, npadx, step, pad_size;

    /**
     * The bin centers, in radians. Internally, they always range from -pi/2 to
     * pi/2.
     */
    protected double[] bins;

    /**
     * The directionality histogram, one array per processor (3 in the case of a
     * ColorProcessor).
     */
    protected ArrayList< double[] > histograms;

    private FloatProcessor padded_square_block;

    private float[] window_pixels;

    /** Store fit results when fit method is called. */
    protected ArrayList< double[] > params_from_fit;

    /** Store goodness of fit results when fit method is called. */
    protected double[] goodness_of_fit;

    /** Store a String representing the fitting function. */
    protected String fit_string;

    /** Used to pass the slice we are currently analyzing. */
    private int slice_index;

    /** This stack stores the orientation map. */
    ImageStack orientation_map;

    /*
     * PLUGIN METHODS
     */

    /**
     * Called when this plugin is launched from ImageJ. This method
     * <ol>
     * <li>grabs the current ImagePlus
     * <li>displays the user dialog and sets setting fields accordingly
     * <li>calls the {@link #computeHistograms()} method, which computes the
     * histograms
     * <li>calls the {@link #fitHistograms()} method, which fits the histograms
     * <li>display the results
     * </ol>
     * <p>
     * If the method is called with String arguments, fields are set according
     * to it, and no dialog are displayed (macro recordable).
     *
     */
    @Override
    public void run()
    {
        //the string argument, for instance "nbins=90, start=-90, end=90, method=gradient"
        String arg = Macro.getOptions();
        // Test if we get an image
        imp = WindowManager.getCurrentImage();
        if ( null == imp )
        {
            IJ.error( "Directionality", "No images are open." );
            return;
        }

        final Roi roi = imp.getRoi();
        if ( null != roi )
            imp = new Duplicator().run( imp, 1, imp.getNSlices() );

        // Non-interactive mode?
        if ( null != arg && arg.length() > 0 )
        {

            // Parse possible macro inputs
            boolean endIsSpecified = false;
            String str = parseArgumentString( arg, "nbins=" );
            if ( null != str )
            {
                try
                {
                    nbins = Integer.parseInt( str );
                }
                catch ( final NumberFormatException nfe )
                {
                    IJ.error( "Directionality: bad argument for number of bins: " + str );
                    return;
                }
            }
            str = parseArgumentString( arg, "start=" );
            if ( null != str )
            {
                try
                {
                    bin_start = Double.parseDouble( str );
                }
                catch ( final NumberFormatException nfe )
                {
                    IJ.error( "Directionality: bad argument for start point: " + str );
                    return;
                }
            }
            str = parseArgumentString( arg, "end=" );
            if ( null != str )
            {
                try
                {
                    bin_end = Double.parseDouble( str );
                    endIsSpecified = true;
                }
                catch ( final NumberFormatException nfe )
                {
                    IJ.error( "Directionality: bad argument for end point: " + str );
                    return;
                }
            }
            str = parseArgumentString( arg, "method=" );
            if ( null != str )
            {
                for ( final AnalysisMethod m : AnalysisMethod.values() )
                {
                    if ( m.toCommandName().equalsIgnoreCase( str ) )
                    {
                        method = m;
                    }
                }
            }
            // if not indicated, the angles cover the full range
            if ( ! endIsSpecified )
            {
                bin_end = bin_start + 180;
            }
        }

        for(AnalysisMethod m :AnalysisMethod.values()){
            if(string_method.equals(m.toString()))
                method = m;
        }

        bin_end = Math.min( bin_end , bin_start + 180 );

        // Launch analysis, this will set the directionality field
        computeHistograms();

        // Fit histograms
        fitHistograms();

        // Display results
        final JFrame plot_frame = plotResults();
        final JFrame data_frame = displayFitAnalysis();

//		int x = Math.max(0, imp.getWindow().getLocation().x - plot_frame.getSize().width);
//		int y = imp.getWindow().getLocation().y;
//		plot_frame.setLocation(x, y);

//		y += plot_frame.getHeight();
//		if (y>Toolkit.getDefaultToolkit().getScreenSize().getHeight()) {
//			y = (int) (0.9 * Toolkit.getDefaultToolkit().getScreenSize().getHeight());
//		}
//		data_frame.setLocation(x, y);
        plot_frame.setLocationRelativeTo( imp.getWindow() );
        data_frame.setLocationRelativeTo( plot_frame );
        plot_frame.setVisible( true );
        data_frame.setVisible( true );

        if ( display_table )
        {
            final ResultsTable table = displayResultsTable();
            table.show( "Directionality histograms for " + imp.getShortTitle() + " (using " + method.toString() + ")" );
        }

        if ( build_orientation_map )
        {
            final ImagePlus imp_map = new ImagePlus( "Orientation map for " + imp.getShortTitle(), orientation_map );
            imp_map.show();
            final ImageCanvas canvas_map = imp_map.getCanvas();
            addColorMouseListener( canvas_map, bin_start, bin_end );
        }

        if ( display_color_wheel )
        {
            final ImagePlus cw = generateColorWheel( bin_start, bin_end );
            cw.show();
            final ImageCanvas canvas_cw = cw.getCanvas();
            addColorMouseListener( canvas_cw, bin_start, bin_end );
        }
    }

    /*
     * PUBLIC METHODS
     */

    /**
     * This method runs the analysis on all slices, and store resulting
     * histograms in the histogram fields. Calling this method resets the
     * aforementioned field.
     */
    public void computeHistograms()
    {
        if ( null == imp )
            return;

        // Reset analysis fields
        params_from_fit = null;
        goodness_of_fit = null;

        // Prepare helper fields
        bins = prepareBins( nbins, bin_start, bin_end );
        switch ( method )
        {
            case FOURIER_COMPONENTS:
                initFourierFields();
                break;
            case LOCAL_GRADIENT_ORIENTATION:
                break;
        }

        // Prepare result holder
        final int n_slices = imp.getStackSize();
        histograms = new ArrayList<>(n_slices * imp.getNChannels());
        if ( build_orientation_map )
        {
            orientation_map = new ImageStack( imp.getWidth(), imp.getHeight() );
        }

        // Loop over each slice
        ImageProcessor ip;
        double[] dir = null;
        for ( int i = 0; i < n_slices; i++ )
        {
            slice_index = i;
            ip = imp.getStack().getProcessor( i + 1 );
            for ( int channel_number = 0; channel_number < ip.getNChannels(); channel_number++ )
            {

                // Convert to float processor
                fip = ip.toFloat( channel_number, fip );

                // Dispatch to specialized method
                switch ( method )
                {
                    case FOURIER_COMPONENTS:
                        dir = fourier_component( fip );
                        break;
                    case LOCAL_GRADIENT_ORIENTATION:
                        dir = local_gradient_orientation( fip );
                        break;
                }

                // Normalize directionality
                double sum = dir[ 0 ];
                for ( int j = 1; j < dir.length; j++ )
                {
                    sum += dir[ j ];
                }
                for ( int j = 0; j < dir.length; j++ )
                {
                    dir[ j ] = dir[ j ] / sum;
                }

                histograms.add( dir );
            }
        }
    }

    /**
     * This method generates a {@link ResultsTable} containing the histogram
     * data for display in ImageJ. It can be used to export the data to a CSV
     * file.
     *
     * @return the result table, which show() method must be called to become
     *         visible.
     */
    public ResultsTable displayResultsTable()
    {
        if ( null == histograms )
            return null;

        final ResultsTable table = new ResultsTable();
        table.setPrecision( 9 );
        final String[] names = makeNames();
        double[] dir;
        for ( int i = 0; i < bins.length; i++ )
        {
            table.incrementCounter();
            table.addValue( "Direction (°)", Math.toDegrees( bins[ i ] ));
            for ( int j = 0; j < names.length; j++ )
            {
                dir = histograms.get( j );
                table.addValue( names[ j ], dir[ i ] );
                final double val = CurveFitter.f( CurveFitter.GAUSSIAN, params_from_fit.get( j ), bins[ i ] );
                table.addValue( names[ j ] + "-fit", val );
            }
        }

        return table;
    }

    /**
     * Return the result of analyzing the gaussian fit of the peak. Results are
     * returned in the shape of an ArrayList of double[], one element per slice.
     * The content of the double arrays is as follow:
     * <ol start=0>
     * <li>gaussian peak center
     * <li>gaussian standard deviation
     * <li>amount, that is: the sum of the histogram data from the gaussian
     * center until SIGMA_NUMBER times its standard deviation away
     * <li>the goodness of fit
     * </ol>
     * The periodic nature of the data is taken into account. For the amount
     * value, the actual values of the histogram are summed, not the values from
     * the fit.
     *
     *
     * @return the fit analysis
     */
    public ArrayList< double[] > getFitAnalysis()
    {
        if ( null == histograms )
            return null;

        final ArrayList< double[] > fit_analysis = new ArrayList<>(histograms.size());
        final double[] gof = getGoodnessOfFit();
        double[] params;
        double[] dir;
        double[] analysis;
        double amount, center, std, xn;

        for ( int i = 0; i < histograms.size(); i++ )
        {
            params = params_from_fit.get( i );
            dir = histograms.get( i );
            analysis = new double[ 4 ];

            amount = 0; // we sum under +/- N*sigma, taking periodicity into
            // account
            center = params[ 2 ];
            std = params[ 3 ];
            for ( int j = 0; j < dir.length; j++ )
            {
                xn = bins[ j ];
                if ( Math.abs( xn - center ) > 90.0 )
                { // too far, we want to shift then
                    if ( xn > center )
                    {
                        xn = xn - 180.0;
                    }
                    else
                    {
                        xn = xn + 180.0;
                    }
                }
                if ( ( xn < ( center - SIGMA_NUMBER * std ) ) || ( xn > ( center + SIGMA_NUMBER * std ) ) )
                {
                    continue;
                }
                amount += dir[ j ];
            }

            analysis[ 0 ] = center;
            analysis[ 1 ] = std;
            analysis[ 2 ] = amount;
            analysis[ 3 ] = gof[ i ];
            fit_analysis.add( analysis );
        }
        return fit_analysis;
    }

    /**
     * This method is called to draw the histograms resulting from image
     * analysis. It reads the result in the histogram list field, and use the
     * JFreeChart library to draw a nice plot window. If the
     * {@code fitResults()} method was called before, the fits are also drawn.
     *
     * @return a {@link JFrame} containing the histogram plots, which
     *         setVisible(boolean) method must be called in order to be
     *         displayed
     */
    public JFrame plotResults()
    {
        final XYSeriesCollection histogram_plots = new XYSeriesCollection();
        final LookupPaintScale lut = createLUT( histograms.size() );
        final String[] names = makeNames();
        XYSeries series;

        final double[] degrees_bins = new double[ nbins ];
        for ( int i = 0; i < degrees_bins.length; i++ )
        {
            degrees_bins[ i ] = Math.toDegrees( bins[ i ] );
        }

        // This is where we shift histograms
        double[] dir;
        for ( int i = 0; i < histograms.size(); i++ )
        {
            dir = histograms.get( i );
            series = new XYSeries( names[ i ] );
            for ( int j = 0; j < nbins; j++ )
            {
                series.add( degrees_bins[ j ], dir[ j ] );
            }
            histogram_plots.addSeries( series );
        }
        histogram_plots.setIntervalWidth( Math.toDegrees( bins[ 1 ] - bins[ 0 ] ) );

        // Create chart with histograms
        final JFreeChart chart = ChartFactory.createHistogram(
                "Directionality histograms",
                "Direction (°)",
                "Amount",
                histogram_plots,
                PlotOrientation.VERTICAL,
                true,
                true,
                false );

        // Set the look of histograms
        final XYPlot plot = ( XYPlot ) chart.getPlot();
        final ClusteredXYBarRenderer renderer = new ClusteredXYBarRenderer( 0.3, false );
        float color_index;
        for ( int i = 0; i < histograms.size(); i++ )
        {
            color_index = ( float ) i / ( float ) ( histograms.size() - 1 );
            renderer.setSeriesPaint( i, lut.getPaint( color_index ) );
        }
        plot.setRenderer( 0, renderer );

        // Draw fit results
        if ( null != params_from_fit )
        {
            // Make new X
            final double[] X = new double[ bins.length * 10 ]; // oversample 10
            // times
            for ( int i = 0; i < X.length; i++ )
            {
                X[ i ] = ( degrees_bins[ 0 ] + ( degrees_bins[ nbins - 1 ] - degrees_bins[ 0 ] ) / X.length * i );
            }
            // Create dataset
            final XYSeriesCollection fits = new XYSeriesCollection();
            XYSeries fit_series;
            double val, xn;
            double[] params;
            for ( int i = 0; i < histograms.size(); i++ )
            { // we have to deal with periodic issue here too
                params = params_from_fit.get( i ).clone();
                fit_series = new XYSeries( names[ i ] );
                for (double x : X) {
                    xn = Math.toRadians(x); // back to radians, for the
                    // fit
                    val = CurveFitter.f(CurveFitter.GAUSSIAN, params, xn);
                    fit_series.add(x, val);
                }
                fits.addSeries( fit_series );
            }
            plot.setDataset( 1, fits );
            plot.setRenderer( 1, new XYLineAndShapeRenderer( true, false ) );
            for ( int i = 0; i < histograms.size(); i++ )
            {
                color_index = ( float ) i / ( float ) ( histograms.size() - 1 );
                plot.getRenderer( 1 ).setSeriesPaint( i, lut.getPaint( color_index ) );
            }

        }
        plot.getDomainAxis().setRange( degrees_bins[ 0 ], degrees_bins[ nbins - 1 ] );

        final ChartPanel chartPanel = new ChartPanel( chart );
        chartPanel.setPreferredSize( new java.awt.Dimension( 500, 270 ) );
        final JFrame window = new JFrame( "Directionality for " + imp.getShortTitle() + " (using " + method.toString() + ")" );
        window.add( chartPanel );
        window.validate();
        window.setSize( new java.awt.Dimension( 500, 270 ) );
        return window;
    }

    /**
     * This method tries to fit a gaussian to the highest peak of each
     * directionality histogram, and store fit results in the
     * <code>params_from_fit</code> field. The goodness of fit will be stored in
     * <code>goodness_of_fit</code>.
     */
    public void fitHistograms()
    {
        if ( null == histograms || histograms.size() == 0)
            return;

        params_from_fit = new ArrayList<>(histograms.size());
        goodness_of_fit = new double[ histograms.size() ];

        // Prepare fitter and function
        CurveFitter fitter = null;

        // Loop over slices
        for ( int i = 0; i < histograms.size(); i++ )
        {
            double[] dir = histograms.get( i );

            fitter = new CurveFitter( bins, dir );

            // Do fit
            fitter.doFit( CurveFitter.GAUSSIAN );
            double[] params = fitter.getParams();
            goodness_of_fit[ i ] = fitter.getFitGoodness();
            params[ 3 ] = Math.abs( params[ 3 ] ); // std is positive
            params_from_fit.add( params );
        }

        fit_string = fitter.getFormula();
    }

    /**
     * This method displays the fit analysis results in a {@link JTable}.
     *
     * @return a {@link JFrame} containing the table; its setVisible(boolean)
     *         method must be called in order to be displayed
     */
    public JFrame displayFitAnalysis()
    {
        if ( null == params_from_fit ) { return null; }
        // Display result
        final String[] column_names = {
                "Slice",
                "Direction (°)",
                "Dispersion (°)",
                "Amount",
                "Goodness" };
        final Object[][] table_data = new Object[ params_from_fit.size() ][ column_names.length ];
        final String[] names = makeNames();
        final ArrayList< double[] > fit_analysis = getFitAnalysis();

        for ( int i = 0; i < table_data.length; i++ )
        {
            double[] analysis = fit_analysis.get( i );
            table_data[ i ][ 0 ] = names[ i ];
            table_data[ i ][ 1 ] = String.format( "%.2f", Math.toDegrees( analysis[ 0 ] ) ); // peak
            // center
            table_data[ i ][ 2 ] = String.format( "%.2f", Math.toDegrees( analysis[ 1 ] ) ); // standard
            // deviation
            table_data[ i ][ 3 ] = String.format( "%.2f", analysis[ 2 ] ); // amount
            table_data[ i ][ 4 ] = String.format( "%.2f", analysis[ 3 ] ); // goodness
            // of
            // fit
        }
        final JTable table = new JTable( table_data, column_names );
        table.setPreferredScrollableViewportSize( new Dimension( 500, 70 ) );

        final JScrollPane scrollPane = new JScrollPane( table );
        table.setAutoResizeMode( JTable.AUTO_RESIZE_ALL_COLUMNS );

        final JPanel table_panel = new JPanel( new GridLayout() );
        table_panel.add( scrollPane );
        final JFrame frame = new JFrame( "Directionality analysis for " + imp.getShortTitle() + " (using " + method.toString() + ")" );

        frame.setDefaultCloseOperation( JFrame.DISPOSE_ON_CLOSE );
        // Create and set up the content pane.
        frame.setContentPane( table_panel );

        // Display the window.
        frame.pack();
        return frame;
    }

    /*
     * SETTERS AND GETTERS
     */

    /**
     * Sets the image for analysis. Calling this method resets the field
     * <code>histograms</code> to null.
     *
     * @param imp
     *            the image.
     */
    public void setImagePlus( final ImagePlus imp )
    {
        this.imp = imp;
        histograms = null;
    }

    /**
     * Gets the image analyzed.
     *
     * @return the image.
     */
    public ImagePlus getImagePlus()
    {
        return imp;
    }

    /**
     * Returns the parameters of the gaussian fit of the main peak in histogram.
     * Results are arranged in an ArrayList of double[], one array per slice
     * analyzed. If the fit was not done prior to this method call, it is
     * called. If the method {@link #computeHistograms()} was not called, null
     * is returned.
     * <p>
     * The double array is organized as follow, for the fitting model y = a +
     * (b-a)*exp(-(x-c)*(x-c)/(2*d*d))
     * <ol start=0>
     * <li>a
     * <li>b
     * <li>x
     * <li>d
     * </ol>
     *
     * @see #getGoodnessOfFit()
     * @see #getHistograms()
     * @see #getBins()
     * @return the fitting parameters.
     */
    public ArrayList< double[] > getFitParameters()
    {
        if ( null == params_from_fit )
        {
            fitHistograms();
        }
        return params_from_fit;
    }

    /**
     * Returns the goodness of fit for the gaussian fit; 1 is good, 0 is bad.
     * One value per slice. If the fit was not done prior to this method call,
     * it is called. If the method {@link #computeHistograms()} was not called,
     * null is returned.
     *
     * @see #getFitParameters()
     * @see #getHistograms()
     * @see #getBins()
     * @return the goodness of fit.
     */
    public double[] getGoodnessOfFit()
    {
        if ( null == params_from_fit )
        {
            fitHistograms();
        }
        return goodness_of_fit;
    }

    /**
     * Returns the directionality histograms as an ArrayList of double[], one
     * array per slice.
     *
     * @see #getBins()
     * @see #getFitParameters()
     * @see #getGoodnessOfFit()
     * @return the directionality histograms; is null if the method
     *         {@link #computeHistograms()} was not called before.
     */
    public ArrayList< double[] > getHistograms()
    {
        return histograms;
    }

    /**
     * Returns the center of the bins for the directionality histograms. They
     * are in degrees.
     *
     * @see #getHistograms()
     * @see #getFitParameters()
     * @see #getGoodnessOfFit()
     * @return the bin centers, in degrees.
     */
    public double[] getBins()
    {
        final double[] degree_bins = new double[ nbins ];
        for ( int i = 0; i < degree_bins.length; i++ )
            degree_bins[ i ] = 180 * bins[ i ]  / Math.PI;

        return degree_bins;
    }

    /**
     * Sets the desired number of bins. This resets the <code>histograms</code>
     * field to null.
     *
     * @param bins
     *            the number of bins.
     */
    public void setBinNumber( final int bins )
    {
        nbins = bins;
        histograms = null;
    }

    /**
     * Returns the current number of bins for this instance.
     *
     * @return the number of bins.
     */
    public int getBinNumber()
    {
        return nbins;
    }

    /**
     * Sets the desired start for the angle bins, in degrees. This resets the
     * <code>histograms</code> field to null.
     *
     * @param start
     *            the bin start.
     */
    public void setBinStart( double start )
    {
        bin_start = start;
        bin_end = bin_start + 180;
        histograms = null;
    }

    /**
     * Returns the current value for angle bin start, in degrees.
     *
     * @return the bin start.
     */
    public double getBinStart()
    {
        return bin_start;
    }

    /**
     * Sets the desired end for the angle bins, in degrees. This resets the
     * <code>histograms</code> field to null.
     *
     * @param bin_start
     *            the bin start.
     * @param bin_end
     *            the bin end.
     */
    public void setBinRange( final double bin_start, final double bin_end )
    {
        this.bin_start = bin_start;
        this.bin_end = Math.min( bin_end , bin_start + 180 );
        histograms = null;
    }

    /**
     * Returns the current value for angle bin end, in degrees.
     *
     * @return the bin end.
     */
    public double getBinEnd()
    {
        return bin_end;
    }

    /**
     * Sets the desired method for analysis. This resets the
     * <code>histograms</code> field to null.
     *
     * @param method
     *            the analysis method.
     *
     * @see AnalysisMethod
     */
    public void setMethod( final AnalysisMethod method )
    {
        this.method = method;
        histograms = null;
    }

    /**
     * Returns the analysis method used by this instance.
     *
     * @return the analysis method.
     */
    public AnalysisMethod getMethod()
    {
        return method;
    }

    /**
     * Sets the debug flag.
     *
     * @param flag
     *            the debug flag.
     */
    public void setDebugFlag( final boolean flag )
    {
        debug = flag;
    }

    /**
     * Sets the build orientation map flag
     *
     * @param flag
     *            whether the orientation map shall be built.
     */
    public void setBuildOrientationMapFlag( final boolean flag )
    {
        build_orientation_map = flag;
    }

    /**
     * Returns the orientation map as an {@link ImageStack}, one slice per slice
     * in the source image. Return null if the orientation map flag was not set,
     * or if computation was not done.
     *
     * @return the orientation map.
     *
     */
    public ImageStack getOrientationMap()
    {
        return orientation_map;
    }

    public void setSobel(int size) { sobel_size = size; }
    public int getSobel() {return sobel_size; }
    /*
     * PRIVATE METHODS
     */


    /**
     * This method is used to initialize variables required for the Fourier
     * analysis after parameters have been set by the user.
     */
    private void initFourierFields()
    {
        if ( null == imp )
            return;

        // Compute dimensions
        width = imp.getWidth();
        height = imp.getHeight();
        if ( width == height )
        {
            npadx = 1;
            npady = 1;
            long_side = width;
            small_side = width;
            step = 0;
        }
        else
        {
            small_side = Math.min( width, height );
            long_side = Math.max( width, height );
            final int npad = long_side / small_side + 1;
            if ( width == long_side )
            {
                npadx = npad;
                npady = 1;
            }
            else
            {
                npadx = 1;
                npady = npad;
            }
            final int delta = ( long_side - small_side );
            step = delta / ( npad - 1 );
        }

        // Computes power of 2 image dimension
        pad_size = 2;
        while ( pad_size < small_side )
            pad_size *= 2;
        padded_square_block = new FloatProcessor( pad_size, pad_size );

        // Prepare windowing
        window = getBlackmanProcessor( small_side, small_side );
        window_pixels = ( float[] ) window.getPixels();

        // Prepare polar coordinates
        r = makeRMatrix( pad_size, pad_size );
        theta = makeThetaMatrix( pad_size, pad_size );

        // Prepare filters
        filters = makeFftFilters();

        if ( debug )
        {
            new ImagePlus( "Angular filters", filters ).show();
        }
    }

    /**
     * This method implements the local gradient orientation analysis method.
     * <p>
     * The gradient is calculated using a 5x5 Sobel filter. The gradient
     * orientation within the given angular range is calculated and put in the histogram.
     * The histogram get the square of the norm as value, so that only strong
     * gradient contribute to it. We use the square of the norm, so that the
     * histogram calculated with this method has the same dimension that with
     * the Fourier method.
     *
     * @see #fourier_component(FloatProcessor)
     *
     */
    private double[] local_gradient_orientation( final FloatProcessor ip )
    {
//		double[] dir = new double[nbins]; // histo with #bins
        final double[] norm_dir = new double[ nbins ]; // histo;
        final FloatProcessor grad_x = ( FloatProcessor ) ip.duplicate();
        final FloatProcessor grad_y = ( FloatProcessor ) ip.duplicate();
        final Convolver convolver = new Convolver();
        final float[][] kernels = getSobel(sobel_size);


        convolver.convolveFloat( grad_x, kernels[0], 5, 5 );
        convolver.convolveFloat( grad_y, kernels[1], 5, 5 );

        final float[] pixels_gx = ( float[] ) grad_x.getPixels();
        final float[] pixels_gy = ( float[] ) grad_y.getPixels();
        final float[] pixels_theta = new float[ pixels_gx.length ];
        final float[] pixels_r = new float[ pixels_gx.length ];

        double norm, max_norm = 0.0;
        double angle, angle_degs, angle_ref;
        double wrapped_start, wrapped_end;
        double range1_min, range1_max, range2_min, range2_max;
        int histo_index;
        float dx, dy;

        // generate range of allowed angles within -90 and +90 degs
        wrapped_start = ( ( bin_start + 90 ) % 180 + 180 ) % 180 - 90;
        wrapped_end = ( ( bin_end + 90 ) % 180 + 180 ) % 180 - 90;

        if ( wrapped_end <= wrapped_start ) {
            range1_min = -90;
            range1_max = wrapped_end;
            range2_min = wrapped_start;
            range2_max = 90;
        }
        else
        {
            range1_min = wrapped_start;
            range1_max = wrapped_end;
            range2_min = wrapped_start;
            range2_max = wrapped_end;
        }


        for ( int i = 0; i < pixels_gx.length; i++ )
        {
            dx = pixels_gx[ i ];
            dy = -pixels_gy[ i ]; // upright orientation
            norm = dx * dx + dy * dy; // We keep the square so as to have the
            // same dimension that Fourier
            // components analysis
            angle = Math.atan( dy / dx );
            angle_degs = angle * 180.0 / Math.PI;

            if ( ( angle_degs >= range1_min && angle_degs <= range1_max )
                    || ( angle_degs >= range2_min && angle_degs <= range2_max ) ) {

                if ( norm > max_norm )
                {
                    max_norm = norm;
                }

                angle_ref = angle * 180.0 / Math.PI;
                while ( angle_ref > bin_end ) {
                    angle_ref -= 180;
                }
                while ( angle_ref < bin_start ) {
                    angle_ref += 180;
                }

                pixels_theta[ i ] = ( float ) ( angle_ref ); // deg,

                pixels_r[ i ] = ( float ) norm;

                histo_index = ( int ) ( ( angle_ref - bin_start ) / ( ( bins[ 1 ] - bins[ 0 ] ) * 180.0 / Math.PI ) );
                if ( histo_index == nbins )
                {
                    histo_index = 0; // circular shift in case of exact vertical
                    // orientation
                }
                norm_dir[ histo_index ] += norm; // we put the norm, the stronger
                // the better
            }
            else
            {
                pixels_theta[ i ] = ( float ) ( bins[ 0 ] * 180. / Math.PI ); // arbitrary value, since norm will be zero
            }
        }

        if ( build_orientation_map )
        {
            final float[] pixels = ( float[] ) ip.getPixels();
            float max_brightness = Float.NEGATIVE_INFINITY;
            float min_brightness = Float.POSITIVE_INFINITY;


            for (float pixel : pixels) {
                if (pixel > max_brightness) {
                    max_brightness = pixel;
                }
                if (pixel < min_brightness) {
                    min_brightness = pixel;
                }
            }
            final ColorProcessor cp = new ColorProcessor( ip.getWidth(), ip.getHeight() );
            final byte[] H = new byte[ pixels_r.length ];
            final byte[] S = new byte[ pixels_r.length ];
            for ( int i = 0; i < pixels_r.length; i++ )
            {
                H[ i ] = ( byte ) ( 255.0 * ( pixels_theta[ i ] - bins[ 0 ] * 180. / Math.PI) / ( bin_end - bin_start )  );
                S[ i ] = ( byte ) ( 255.0 * pixels_r[ i ] / max_norm ); // Math.log10(1.0
                // +
                // 9.0*pixels_r[i]
                // /
                // max_norm)
                // );
            }
            final byte[] B = ( byte[] ) ip.convertToByte( true ).getPixels();
            cp.setHSB( H, S, B );
            orientation_map.addSlice( makeNames()[ slice_index ], cp );
        }

        return norm_dir;
    }

    /**
     * This method implements the Fourier component analysis method. The method
     * {@link #initFourierFields()} must be called before this method is.
     * <p>
     * Images are chopped in squares, and the Fourier spectrum of each square is
     * calculated. For the spectrum, we get the angular component using the
     * filters generated by {@link #makeFftFilters()}.
     * </p>
     * <p>
     * We return the results as a double array, containing the amount of
     * orientation for each angles specified in the {@link #bins} field.
     * </p>
     * <p>
     * See {@link #local_gradient_orientation(FloatProcessor)}
     * </p>
     */
    private double[] fourier_component( final FloatProcessor ip )
    {
        final Roi original_square = new Roi( ( pad_size - small_side ) / 2, ( pad_size - small_side ) / 2, small_side, small_side );

        float[] fpx, spectrum_px;
        final double[] dir = new double[ nbins ];
        ImageProcessor square_block;
        Roi square_roi;
        FHT fft, pspectrum;
        FloatProcessor small_pspectrum;

        ImageStack spectra = null;
        if ( debug )
        {
            spectra = new ImageStack( small_side, small_side );
        }

        FloatProcessor[] hue_arrays = null, saturation_arrays = null;
        if ( build_orientation_map )
        {
            hue_arrays = new FloatProcessor[ npadx * npady ];
            saturation_arrays = new FloatProcessor[ npadx * npady ];
        }

        // Overall maximum of the weights
        float max_norm = 0.0f;
        // If the image is not square, split it in small square padding all the
        // image
        for ( int ix = 0; ix < npadx; ix++ )
        {

            for ( int iy = 0; iy < npady; iy++ )
            {

                // Extract a square block from the image
                square_roi = new Roi( ix * step, iy * step, small_side, small_side );
                ip.setRoi( square_roi );
                square_block = ip.crop();

                // Window the block
                final float[] block_pixels = ( float[] ) square_block.getPixels();
                for ( int i = 0; i < block_pixels.length; i++ )
                {
                    block_pixels[ i ] *= window_pixels[ i ];
                }

                // Pad the block with a power of 2 size
                padded_square_block.setValue( 0.0 );
                padded_square_block.fill();
                padded_square_block.insert( square_block, ( pad_size - small_side ) / 2, ( pad_size - small_side ) / 2 );

                // Computes its FFT
                fft = new FHT( padded_square_block );
                fft.setShowProgress( false );
                fft.transform();
                fft.swapQuadrants();

                // Get a centered power spectrum with right size
                pspectrum = fft.conjugateMultiply( fft );
                pspectrum.setRoi( original_square );
                spectrum_px = ( float[] ) pspectrum.getPixels(); // small_pspectrum.getPixels();

                if ( debug )
                {
                    assert spectra != null;
                    small_pspectrum = ( FloatProcessor ) pspectrum.crop();
                    spectra.addSlice( "block nbr " + ( ix + 1 ) * ( iy + 1 ), displayLog( small_pspectrum ) );
                }

                // For orientation map
                float[] weights = null, max_weights = null, best_angle = null;
                FHT tmp;
                FloatProcessor small_tmp;
                float[] tmp_px, small_tmp_px;

                if ( build_orientation_map )
                {
                    weights = new float[ small_side * small_side ];
                    max_weights = new float[ small_side * small_side ];
                    best_angle = new float[ small_side * small_side ];
                }

                // Loop over all bins
                for ( int bin = 0; bin < nbins; bin++ )
                {

                    // Get filter pixels
                    fpx = ( float[] ) filters.getPixels( bin + 1 );

                    // Loop over all pixels
                    if ( build_orientation_map )
                    {

                        tmp = fft.getCopy();
                        tmp.setShowProgress( false );
                        tmp_px = ( float[] ) tmp.getPixels();
                        for ( int i = 0; i < spectrum_px.length; i++ )
                        {
                            // Computes angular density
                            dir[ bin ] += spectrum_px[ i ] * fpx[ i ]; // will
                            // sum
                            // out
                            // with
                            // every
                            // block
                            // Build orientation map if needed
                            tmp_px[ i ] *= fpx[ i ];
                        }
                        tmp.inverseTransform();
                        tmp.setRoi( original_square );
                        small_tmp = ( FloatProcessor ) tmp.crop();

                        // Build angular statistics arrays -> 2nd loop
                        small_tmp_px = ( float[] ) small_tmp.getPixels();
                        assert weights != null; // Asserstion w.r.t. regression
                        for ( int j = 0; j < small_tmp_px.length; j++ )
                        {

                            weights[ j ] = small_tmp_px[ j ] * small_tmp_px[ j ];
                            if ( weights[ j ] > max_weights[ j ] )
                            {
                                max_weights[ j ] = weights[ j ];
                                best_angle[ j ] = ( float ) bins[ bin ]; // radians
                            }
                            // Overall maximum calculation
                            if ( weights[ j ] > max_norm )
                            {
                                max_norm = weights[ j ];
                            }
                        }

                    }
                    else
                    {

                        for ( int i = 0; i < spectrum_px.length; i++ )
                        {
                            // Computes angular density, and that's all
                            dir[ bin ] += spectrum_px[ i ] * fpx[ i ]; // will
                            // sum
                            // out
                            // with
                            // every
                            // block
                        }
                    }

                } // end loop over all bins

                // Store results
                if ( build_orientation_map )
                {
                    assert hue_arrays != null; // Sanity checking w.r.t. regression
                    assert saturation_arrays != null;
                    hue_arrays[ ix + npadx * iy ] = new FloatProcessor( ip.getWidth(), ip.getHeight() );
                    hue_arrays[ ix + npadx * iy ].insert( new FloatProcessor( small_side, small_side, best_angle, null ), ix * step, iy * step );
                    saturation_arrays[ ix + npadx * iy ] = new FloatProcessor( ip.getWidth(), ip.getHeight() );
                    saturation_arrays[ ix + npadx * iy ].insert( new FloatProcessor( small_side, small_side, max_weights, null ), ix * step, iy * step );
                }
            }
        }

        // Reconstruct final orientation map
        if ( build_orientation_map )
        {
            final FloatProcessor big_hue = new FloatProcessor( ip.getWidth(), ip.getHeight() );
            final FloatProcessor big_saturation = new FloatProcessor( ip.getWidth(), ip.getHeight() );
            final float[] big_hue_px = ( float[] ) big_hue.getPixels();
            final float[] big_saturation_px = ( float[] ) big_saturation.getPixels();
            float[] saturation_px, hue_px;
            for ( int ix = 0; ix < npadx; ix++ )
            {
                for ( int iy = 0; iy < npady; iy++ )
                {
                    hue_px = ( float[] ) hue_arrays[ ix + npadx * iy ].getPixels();
                    saturation_px = ( float[] ) saturation_arrays[ ix + npadx * iy ].getPixels();
                    for ( int i = 0; i < big_hue_px.length; i++ )
                    {
                        if ( ( 255 * saturation_px[ i ] / max_norm ) >= big_saturation_px[ i ] )
                        {
                            big_saturation_px[ i ] = ( 255 * saturation_px[ i ] / max_norm );
                            big_hue_px[ i ] = ( float ) ( 255 * ( ( ( hue_px[ i ] - bins[ 0 ] ) / ( bins[ nbins - 1 ] - bins[ 0 ] ) ) ) );
                        }
                    }
                }
            }

            final ByteProcessor big_brightness = ( ByteProcessor ) ip.convertToByte( true );
            final ColorProcessor cp = new ColorProcessor( ip.getWidth(), ip.getHeight() );
            cp.setHSB(
                    ( byte[] ) big_hue.convertToByte( false ).getPixels(),
                    ( byte[] ) big_saturation.convertToByte( false ).getPixels(),
                    ( byte[] ) big_brightness.getPixels() );
            orientation_map.addSlice( makeNames()[ slice_index ], cp );
        }

        if ( debug )
        {
            new ImagePlus( "Log10 power FFT of " + makeNames()[ slice_index ], spectra ).show();
        }

        return dir;
    }

    /**
     * This method generates the angular filters used by the Fourier analysis.
     * It reads the fields {@link #nbins}, {@link #bin_start} to determine how
     * many individual angle filter to generate, and {@link #pad_size} to
     * determine the image filter size. As such, they must be set before calling
     * this method.
     *
     * <p>
     * See {@link #fourier_component(FloatProcessor)}, {@link #prepareBins(int, double, double)}
     * </p>
     *
     * @return an {@link ImageStack} made of each individual angular filter
     */
    private ImageStack makeFftFilters()
    {
        final ImageStack filters = new ImageStack( pad_size, pad_size, nbins );
        float[] pixels;

        final float[] r_px = ( float[] ) r.getPixels();
        final float[] theta_px = ( float[] ) theta.getPixels();

        double current_r, current_theta, theta_c, angular_part, radial_part;
        final double theta_bw = bins[ 1 ] - bins[ 0 ];
        final double r_c = pad_size / 4.;
        final double r_bw = r_c / 2;

        for ( int i = 1; i <= nbins; i++ )
        {

            pixels = new float[ pad_size * pad_size ];

            theta_c = bins[ i - 1];
            while ( theta_c >= Math.PI ) {
                theta_c -= Math.PI;
            }
            while ( theta_c < 0 ) {
                theta_c += Math.PI;
            }
            theta_c -= Math.PI/2;

            for ( int index = 0; index < pixels.length; index++ )
            {

                current_r = r_px[ index ];
                if ( current_r < FREQ_THRESHOLD || current_r > pad_size / 2. )
                {
                    continue;
                }
                radial_part = Math.exp( -( current_r - r_c ) * ( current_r - r_c ) / ( r_bw * r_bw ) );

                current_theta = theta_px[ index ];
                if ( Math.abs( current_theta - theta_c ) < theta_bw )
                {
                    angular_part = Math.cos( ( current_theta - theta_c ) / theta_bw * Math.PI / 2.0 );
                    angular_part = angular_part * angular_part;

                }
                else if ( Math.abs( current_theta - ( theta_c - Math.PI ) ) < theta_bw )
                {
                    angular_part = Math.cos( ( current_theta - ( theta_c - Math.PI ) ) / theta_bw * Math.PI / 2.0 );
                    angular_part = angular_part * angular_part;
                }
                else if ( Math.abs( current_theta - ( theta_c + Math.PI ) ) < theta_bw )
                {
                    angular_part = Math.cos( ( current_theta - ( theta_c + Math.PI ) ) / theta_bw * Math.PI / 2.0 );
                    angular_part = angular_part * angular_part;
                }
                else
                {
                    continue; // leave it to 0
                }

                pixels[ index ] = ( float ) ( angular_part * radial_part );

            }

            filters.setPixels( pixels, i );
            filters.setSliceLabel( "Angle: " + String.format( "%.1f", bins[ i - 1] * 180 / Math.PI ), i );
        }
        return filters;
    }

    /**
     * This method generate a name for each analyzed slice, to display in result
     * tables.
     *
     * @return a String array with the names
     */
    private String[] makeNames()
    {
        final int n_slices = imp.getStack().getSize();
        String[] names;
        String label;
        if ( imp.getType() == ImagePlus.COLOR_RGB )
        {
            names = new String[ 3 * n_slices ];
            for ( int i = 0; i < n_slices; i++ )
            {
                label = imp.getStack().getShortSliceLabel( i + 1 );
                if ( null == label )
                {
                    names[ i * 3] = "Slice_" + ( i + 1 ) + "R";
                    names[ 1 + i * 3 ] = "Slice_" + ( i + 1 ) + "G";
                    names[ 2 + i * 3 ] = "Slice_" + ( i + 1 ) + "B";
                }
                else
                {
                    names[ i * 3] = label + "_R";
                    names[ 1 + i * 3 ] = label + "_G";
                    names[ 2 + i * 3 ] = label + "_B";
                }
            }
        }
        else
        {
            if ( n_slices <= 1 ) { return new String[] { imp.getShortTitle() }; }
            names = new String[ n_slices ];
            for ( int i = 0; i < n_slices; i++ )
            {
                label = imp.getStack().getShortSliceLabel( i + 1 );
                if ( null == label )
                {
                    names[ i ] = "Slice_" + ( i + 1 );
                }
                else
                {
                    names[ i ] = label;
                }
            }
        }
        return names;
    }

    /*
     * STATIC METHODS
     */

    public static ImagePlus generateColorWheel( final double angle_start, final double angle_end )
    {
        final int cw_height = 256;
        final int cw_width = cw_height / 2;
        final int offset = 64;
        final ColorProcessor color_ip = new ColorProcessor( cw_width + offset, cw_height );
        FloatProcessor R = makeRMatrix( cw_height + 2 * offset, cw_height );
        FloatProcessor T = makeThetaMatrix( cw_height + 2 * offset, cw_height );
        final Roi half_roi = new Roi( cw_height / 2 + offset, 0, cw_width + offset, cw_height );
        R.setRoi( half_roi );
        R = ( FloatProcessor ) R.crop();
        T.setRoi( half_roi );
        T = ( FloatProcessor ) T.crop();
        final float[] r = ( float[] ) R.getPixels();
        final float[] t = ( float[] ) T.getPixels();
        final byte[] hue = new byte[ r.length ];
        final byte[] sat = new byte[ r.length ];
        final byte[] bgh = new byte[ r.length ];
        for ( int i = 0; i < t.length; i++ )
        {
            if ( r[ i ] > cw_width )
            {
                hue[ i ] = 0;
                sat[ i ] = 0;
            }
            else
            {
                hue[ i ] = ( byte ) ( 255 * ( t[ i ] + Math.PI / 2 ) / Math.PI );
                sat[ i ] = ( byte ) ( 255 * r[ i ] / cw_width );
            }
            bgh[ i ] = ( byte ) 255;
        }
        color_ip.setHSB( hue, sat, bgh );
        color_ip.setBackgroundValue( 255 );
        color_ip.filterRGB(ColorProcessor.RGB_TRANSLATE, offset);
        color_ip.setColor( Color.WHITE );
        final Roi fill_roi = new Roi( 0, 0, offset, cw_height );
        color_ip.setRoi( fill_roi );
        color_ip.fill( fill_roi );

        color_ip.setColor( Color.BLACK );
        color_ip.setJustification( ImageProcessor.RIGHT_JUSTIFY );
        final String txt_min = String.valueOf( angle_start );
        final String txt_mid = String.valueOf( 0.5 * ( angle_start + angle_end ) );
        final String txt_max = String.valueOf( angle_end );
        final TextRoi text_roi_min = new TextRoi( offset - 5, cw_height - 25, txt_min );
        text_roi_min.drawPixels( color_ip );
        final TextRoi text_roi_max = new TextRoi( offset - 5, 0, txt_max );
        text_roi_max.drawPixels( color_ip );
        final TextRoi text_roi_mid = new TextRoi( offset - 5, cw_height / 2 - 13, txt_mid );
        text_roi_mid.drawPixels( color_ip );

        return new ImagePlus( "Color wheel", color_ip );
    }

    protected static void addColorMouseListener( final ImageCanvas canvas, final double angle_start, final double angle_end )
    {

        final MouseMotionListener ml = new MouseMotionListener()
        {
            @Override
            public void mouseDragged( final MouseEvent e )
            {}

            @Override
            public void mouseMoved( final MouseEvent e )
            {
                final Point coord = canvas.getCursorLoc();
                final int x = coord.x;
                final int y = coord.y;
                try
                {
                    final ColorProcessor cp = ( ColorProcessor ) canvas.getImage().getProcessor();
                    final int c = cp.getPixel( x, y );
                    final int r = ( c & 0xff0000 ) >> 16;
                    final int g = ( c & 0xff00 ) >> 8;
                    final int b = c & 0xff;
                    final float[] hsb = Color.RGBtoHSB( r, g, b, null );
                    final float angle = (float) (hsb[ 0 ] * (angle_end - angle_start) + angle_start);
                    final float amount = hsb[ 1 ];
                    IJ.showStatus( String.format( "Orientation: %5.1f ° - Amont: %5.1f %%", angle, 100 * amount ) );
                }
                catch ( final ClassCastException ignored)
                {
                }
            }
        };
        canvas.addMouseMotionListener( ml );
    }

    /**
     * Generate a bin array of angle in radians, from -pi/2 to pi/2.
     *
     * @return a double array of n elements, the angles in radians
     * @param n
     *            the number of elements to generate
     * @param first
     *            the angle of the first bin, in degrees
     * @param last
     *            the angle of the last bin, in degrees
     */
    protected static double[] prepareBins( final int n, final double first, final double last )
    {
        final double[] bins = new double[ n ];
        for ( int i = 0; i < n; i++ )
        {
            bins[ i ] = ( first + i * ( last - first) / ( n - 1) ) / 180. * Math.PI ;
        }
        return bins;
    }

    /**
     * Utility method to analyze the content of the argument string passed by
     * ImageJ to this plugin using the {@code setup(String, ImagePlus)} method.
     * Not as clever as it could be.
     *
     * @param argument_string
     *            the argument string to parse
     * @param command_str
     *            the command to search for
     * @return a string containing the value after the command string, null if
     *         the command string was not found
     */
    protected static String parseArgumentString( final String argument_string, final String command_str )
    {
        if ( argument_string.contains( command_str ) )
        {
            final int narg = argument_string.indexOf( command_str ) + command_str.length();
            int next_arg = argument_string.indexOf( ",", narg );
            if ( next_arg == -1 )
            {
                next_arg = argument_string.length();
            }
            return argument_string.substring( narg, next_arg );
        }
        return null;
    }

    /**
     * This utility method returns a FloatProcessor with the log10 of each pixel
     * in the {@link FloatProcessor} given in argument. Usefull to display power
     * spectrum.
     *
     * @param ip
     *            the source FloatProcessor
     * @return a new {@link FloatProcessor}
     */
    protected static FloatProcessor displayLog( final FloatProcessor ip )
    {
        final FloatProcessor log10 = new FloatProcessor( ip.getWidth(), ip.getHeight() );
        final float[] log10_pixels = ( float[] ) log10.getPixels();
        final float[] pixels = ( float[] ) ip.getPixels();
        for ( int i = 0; i < pixels.length; i++ )
        {
            log10_pixels[ i ] = ( float ) Math.log10( 1 + pixels[ i ] );
        }
        return log10;
    }

    /**
     * This utility generates a <b>periodic</b> Blackman window over n points.
     *
     * @param n
     *            the number of point in the window
     * @return a double array containing the Blackman window
     * @see #getBlackmanProcessor(int, int)
     */
    protected static double[] getBlackmanPeriodicWindow1D( final int n )
    {
        final double[] window = new double[ n ];
        for ( int i = 0; i < window.length; i++ )
        {
            window[ i ] = 0.42 - 0.5 * Math.cos( 2 * Math.PI * i / n ) + 0.08 * Math.cos( 4 * Math.PI / n );
        }
        return window;
    }

    /**
     * Generate a 2D Blackman window used in this plugin before computing FFT,
     * so as to avoid the cross artifact at x=0 and y=0.
     *
     * @param nx
     *            the width in pixel of the desired window
     * @param ny
     *            the hirght in pixel of the desired window
     * @return the window, as a FloatProcessor
     * @see #getBlackmanPeriodicWindow1D(int)
     */
    protected static FloatProcessor getBlackmanProcessor( final int nx, final int ny )
    {
        final FloatProcessor bpw = new FloatProcessor( nx, ny );
        final float[] pixels = ( float[] ) bpw.getPixels();
        final double[] bpwx = getBlackmanPeriodicWindow1D( nx );
        final double[] bpwy = getBlackmanPeriodicWindow1D( nx );
        int ix, iy;
        for ( int i = 0; i < pixels.length; i++ )
        {
            iy = i / nx;
            ix = i % nx;
            pixels[ i ] = ( float ) ( bpwx[ ix ] * bpwy[ iy ] );
        }
        return bpw;
    }

    /**
     * Generate a 2D matrix of the radius polar coordinates, centered in the
     * middle of the image.
     *
     * @param nx
     *            the width in pixel of the desired matrix
     * @param ny
     *            the height in pixel of the desired matrix
     * @return the coordinate matrix, as a FloatProcessor, where values range
     *         from 0 to sqrt(nx^2+ny^2)/2
     * @see #makeThetaMatrix(int, int)
     */
    protected static FloatProcessor makeRMatrix( final int nx, final int ny )
    {
        final FloatProcessor r = new FloatProcessor( nx, ny );
        final float[] pixels = ( float[] ) r.getPixels();
        final int xc = nx / 2;
        final int yc = ny / 2;
        int ix, iy;
        for ( int i = 0; i < pixels.length; i++ )
        {
            iy = i / nx;
            ix = i % nx;
            pixels[ i ] = ( float ) Math.sqrt( ( ix - xc ) * ( ix - xc ) + ( iy - yc ) * ( iy - yc ) );
        }
        return r;
    }

    /**
     * Generate a 2D matrix of the angle polar coordinates, centered in the
     * middle of the image.
     *
     * @param nx
     *            the width in pixel of the desired matrix
     * @param ny
     *            the height in pixel of the desired matrix
     * @return the coordinate matrix, as a FloatProcessor, where angles are in
     *         radians, and range from -pi to pi
     * @see #makeRMatrix(int, int)
     */
    protected static FloatProcessor makeThetaMatrix(final int nx, final int ny )
    {
        final FloatProcessor theta = new FloatProcessor( nx, ny );
        final float[] pixels = ( float[] ) theta.getPixels();
        final int xc = nx / 2;
        final int yc = ny / 2;
        int ix, iy;
        for ( int i = 0; i < pixels.length; i++ )
        {
            iy = i / nx;
            ix = i % nx;
            pixels[ i ] = ( float ) Math.atan2( -( iy - yc ), ix - xc ); // so
            // that
            // we
            // have
            // upright
            // orientation
        }
        return theta;
    }

    /**
     * Generate a bluish to greenish to redish LUT for the display of
     * histograms.
     *
     * @param ncol
     *            the number of colors in the LUT
     * @return the LUT
     */
    protected static LookupPaintScale createLUT(final int ncol )
    {
        final float[][] colors = new float[][] {
                { 0, 75 / 255f, 150 / 255f },
                { 0.1f, 0.8f, 0.1f },
                { 150 / 255f, 75 / 255f, 0 }
        };
        final float[] limits = new float[] { 0, 0.5f, 1 };
        final LookupPaintScale lut = new LookupPaintScale( 0, 1, Color.BLACK );
        float val;
        float r, g, b;
        for ( int j = 0; j < ncol; j++ )
        {
            val = j / ( ncol - 0.99f );
            int i;
            for ( i = 0; i < limits.length; i++ )
            {
                if ( val < limits[ i ] )
                {
                    break;
                }
            }
            i = i - 1;
            float v = (val - limits[i]) / (limits[i + 1] - limits[i]);
            r = colors[ i ][ 0 ] + v * ( colors[ i + 1 ][ 0 ] - colors[ i ][ 0 ] );
            g = colors[ i ][ 1 ] + v * ( colors[ i + 1 ][ 1 ] - colors[ i ][ 1 ] );
            b = colors[ i ][ 2 ] + v * ( colors[ i + 1 ][ 2 ] - colors[ i ][ 2 ] );
            lut.add( val, new Color( r, g, b ) );
        }
        return lut;
    }

    /**
     * Generate a sobel filter in the x and y direction for the given size
     * Method adapted from https://stackoverflow.com/questions/9567882/sobel-filter-kernel-of-large-size
     *
     * @param size size of the sobel filter as a square
     *
     * @return the two sobel kernels, x and y as floats
     **/
    protected float[][] getSobel(int size) {
        if(size%2 == 0) {
            logService.info("Sobel filter size cannot be even. Adding one to make odd");
            size += 1;
        }
        float[] kernel_y = new float[size*size];
        float[] kernel_x = new float[size*size];
        int[] center = new int[]{size/2, size/2};

        // X kernel
        for(int y = 0; y < size; y++) {
            for(int x = 0; x < size; x++) {
                int index = y * size + x;
                kernel_x[index] = index;
            }
        }

        // X kernel
        for(int y = 0; y < size; y++) {
            for(int x = 0; x < size; x++) {
                int index = y * size + x;
                kernel_y[index] = index;
            }
        }
        kernel_y = new float[] {
                -5f,  -4f, 0f,  4f,  5f,
                -8f, -10f, 0f, 10f,  8f,
                -10f, -20f, 0f, 20f, 10f,
                -8f, -10f, 0f, 10f,  8f,
                -5f,  -4f, 0f,  4f,  5f };
        kernel_x = new float[] {
                5f,   8f,  10f,   8f,  5f,
                4f,  10f,  20f,  10f,  4f,
                0f,   0f,   0f,   0f,  0f,
                -4f, -10f, -20f, -10f, -4f,
                -5f,  -8f, -10f,  -8f, -5f };

        return new float[][] {kernel_x, kernel_y};
    }

    /*
     * MAIN METHOD
     */

    public static void main( final String[] args )
    {

        // Generate a test image
        final ImagePlus imp = NewImage.createShortImage( "Lines", 400, 400, 1, NewImage.FILL_BLACK );
        final ImageProcessor ip = imp.getProcessor();
        ip.setLineWidth( 4 );
        ip.setColor( Color.WHITE );
        // line, 30º
        final Line line_30deg = new Line( 10.0, 412.0, 446, 160.2753 );
        // line, 30º
        final Line line_30deg2 = new Line( 10.0, 312.0, 446.4102, 60.2753 );
        // line, -60º
        final Line line_m60deg = new Line( 10.0, 10, 300.0, 512.2947 );

        final Line[] rois = new Line[] { line_30deg, line_30deg2, line_m60deg };
        for ( final Line roi : rois )
        {
            roi.setStrokeWidth(2);
            ip.draw( roi );
        }
        final GaussianBlur smoother = new GaussianBlur();
        smoother.blurGaussian( ip, 2.0, 2.0, 1e-2 );
        imp.show();

        AnalysisMethod method;
        ArrayList< double[] > fit_results;
        double center;

        final OwnDirectionality da = new OwnDirectionality();
        da.setImagePlus( imp );

        da.setBinNumber( 60 );
        da.setBinStart( -90 );
//		da.setBinRange( -80, 60 );

        da.setBuildOrientationMapFlag( false );
        da.setDebugFlag( false );

        method = AnalysisMethod.FOURIER_COMPONENTS;
//		method = AnalysisMethod.LOCAL_GRADIENT_ORIENTATION;
        da.setMethod( method );
        da.computeHistograms();
        fit_results = da.getFitParameters();
        center = fit_results.get( 0 )[ 2 ];
        System.out.println( "With method: " + method );
        System.out.printf("Found maxima at %.1f, expected it at 30°.\n%n", center/Math.PI*180.);
//		new ImagePlus("Orientation map for "+imp.getShortTitle(),da.getOrientationMap()).show();

        /*
         * method = AnalysisMethod.LOCAL_GRADIENT_ORIENTATION;
         * da.setMethod(method); da.computesHistograms(); fit_results =
         * da.getFitParameters(); center = fit_results.get(0)[2];
         * System.out.println("With method: "+method);
         * System.out.println(String.
         * format("Found maxima at %.1f, expected it at 30º.\n", center, 30));
         * new ImagePlus("Orientation map for "+imp.getShortTitle(),da.
         * getOrientationMap()).show();
         */

//		ImagePlus cw = generateColorWheel( da.getBinStart(), da.getBinEnd() );
//		cw.show();
//		addColorMouseListener(cw.getCanvas(), da.getBinStart(), da.getBinEnd() );

        da.plotResults().setVisible( true );
        da.displayResultsTable().show( "Table" );

    }

}

