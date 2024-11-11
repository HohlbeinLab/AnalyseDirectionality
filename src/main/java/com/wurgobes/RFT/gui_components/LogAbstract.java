package com.wurgobes.RFT.gui_components;

/**
 * This class and all other classes in this folder were adapted from (<a href="https://github.com/Biomedical-Imaging-Group/OrientationJ">OrientationJ</a>)
 * Please refer to their GitHub page for further documentation.
 */
public interface LogAbstract {

    /**
     * Set a value and a message in the progress bar.
     */
    void progress(String msg, int value);

    /**
     * Set a value and a message in the progress bar.
     */
    void increment(double inc);

    /**
     * Set a value in the progress bar.
     */
    void setValue(int value);

    /**
     * Set a message in the progress bar.
     */
    void setMessage(String msg);

    /**
     * Set a value and a message in the progress bar.
     */
    void progress(String msg, double value);
    /**
     * Set to 0 the progress bar.
     */
    void reset();

    /**
     * Set to 100 the progress bar.
     */
    void finish();

    /**
     * Set to 100 the progress bar with an additional message.
     */
    void finish(String msg);

}

