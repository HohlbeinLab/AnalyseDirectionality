package com.wurgobes.RFT.gui_components;

import javax.swing.*;
import java.awt.*;

/**
 * This class and all other classes in this folder were adapted from (<a href="https://github.com/Biomedical-Imaging-Group/OrientationJ">OrientationJ</a>)
 * Please refer to their GitHub page for further documentation.
 * <p>
 * This class extends the JPanel to create grid panel given the possibility to
 * place Java compoments in an organized manner in the dialog box.
 *
 * @author Daniel Sage, Biomedical Imaging Group, EPFL, Lausanne, Switzerland.
 *
 */
public class GridPanel extends JPanel {

    private GridBagLayout		layout			= new GridBagLayout();
    private GridBagConstraints	constraint		= new GridBagConstraints();
    private int					defaultSpace	= 3;

    /**
     * Constructor.
     */
    public GridPanel() {
        super();
        setLayout(layout);
        setBorder(BorderFactory.createEtchedBorder());
    }

    /**
     * Constructor.
     */
    public GridPanel(int defaultSpace) {
        super();
        setLayout(layout);
        this.defaultSpace = defaultSpace;
        setBorder(BorderFactory.createEtchedBorder());
    }

    /**
     * Constructor.
     */
    public GridPanel(boolean border) {
        super();
        setLayout(layout);
        if (border) {
            setBorder(BorderFactory.createEtchedBorder());
        }
    }

    /**
     * Constructor.
     */
    public GridPanel(String title) {
        super();
        setLayout(layout);
        setBorder(BorderFactory.createTitledBorder(title));
    }

    /**
     * Constructor.
     */
    public GridPanel(boolean border, int defaultSpace) {
        super();
        setLayout(layout);
        this.defaultSpace = defaultSpace;
        if (border) {
            setBorder(BorderFactory.createEtchedBorder());
        }
    }

    /**
     * Constructor.
     */
    public GridPanel(String title, int defaultSpace) {
        super();
        setLayout(layout);
        this.defaultSpace = defaultSpace;
        setBorder(BorderFactory.createTitledBorder(title));
    }

    /**
     * Specify the defaultSpace.
     */
    public void setSpace(int defaultSpace) {
        this.defaultSpace = defaultSpace;
    }

    /**
     * Place a component in the northwest of the cell.
     */
    public void place(int row, int col, JComponent comp) {
        place(row, col, 1, 1, defaultSpace, comp);
    }

    /**
     * Place a component in the northwest of the cell.
     */
    public void place(int row, int col, int space, JComponent comp) {
        place(row, col, 1, 1, space, comp);
    }

    /**
     * Place a component in the northwest of the cell.
     */
    public void place(int row, int col, int width, int height, JComponent comp) {
        place(row, col, width, height, defaultSpace, comp);
    }

    /**
     * Place a component in the northwest of the cell.
     */
    public void place(int row, int col, int width, int height, int space, JComponent comp) {
        constraint.gridx = col;
        constraint.gridy = row;
        constraint.gridwidth = width;
        constraint.gridheight = height;
        constraint.anchor = GridBagConstraints.NORTHWEST;
        constraint.insets = new Insets(space, space, space, space);
        constraint.fill = GridBagConstraints.HORIZONTAL;
        layout.setConstraints(comp, constraint);
        add(comp);
    }

}
