package com.wurgobes.RFT.gui_components;

import javax.swing.*;
import javax.swing.text.BadLocationException;
import javax.swing.text.Document;
import java.awt.*;

/**
 * This class and all other classes in this folder were adapted from (<a href="https://github.com/Biomedical-Imaging-Group/OrientationJ">OrientationJ</a>)
 * Please refer to their GitHub page for further documentation.
 * <p>
 * This class extends the Java JEditorPane to make a easy to use panel to
 * display HTML information.
 *
 * @author Daniel Sage, Biomedical Imaging Group, EPFL, Lausanne, Switzerland.
 *
 */
public class Credits extends JEditorPane {

    private String	html		= "";
    private String	header		= "";
    private String	footer		= "";
    private String	font		= "verdana";
    private String	color		= "#222222";
    private String	background	= "#f8f8f8";


    private String ref1 =
            "M.I. Gobes et al., "
                    + "Rotated Fourier transform (RFT) enables the quantification of anisotropic structure in high moisture plant-protein extrudates, "
                    + "In Preparation, Food Structure, 2024.";

    private String ref0 =
            "Z. Püspöki et al., "
                    + "Transforms and Operators for Directional Bioimage Analysis: A Survey, "
                    + "Advances in Anatomy, Embryology and Cell Biology, vol. 219, "
                    + "Springer International Publishing, 2016.";

    private String bib0 =
            "@INCOLLECTION(http://bigwww.epfl.ch/publications/puespoeki1603.html,\n" +
                    "AUTHOR=\"P{\\\"{u}}sp{\\\"{o}}ki, Z. and Storath, M. and Sage, D. and Unser,\n" +
                    "        M.\",\n" +
                    "TITLE=\"Transforms and Operators for Directional Bioimage Analysis: {A}\n" +
                    "        Survey\",\n" +
                    "BOOKTITLE=\"Focus on Bio-Image Informatics\",\n" +
                    "PUBLISHER=\"Springer International Publishing\",\n" +
                    "YEAR=\"2016\",\n" +
                    "editor=\"De Vos, W.H. and Munck, S. and Timmermans, J.-P.\",\n" +
                    "volume=\"219\",\n" +
                    "series=\"Advances in Anatomy, Embryology and Cell Biology\",\n" +
                    "type=\"\",\n" +
                    "chapter=\"3\",\n" +
                    "pages=\"69--93\",\n" +
                    "address=\"\",\n" +
                    "edition=\"\",\n" +
                    "month=\"May 21,\",\n" +
                    "note=\"\")";

    public Credits() {
        header += "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 3.2//EN\">\n";
        header += "<html><head>\n";
        header += "<style>body {background-color:" + background + "; color:" + color + "; font-family: " + font + ";margin:4px}</style>\n";
        header += "<style>h1 {color:#555555; font-size:1.0em; font-weight:bold; padding:1px; margin:1px; padding-top:5px}</style>\n";
        header += "<style>h2 {color:#333333; font-size:0.9em; font-weight:bold; padding:1px; margin:1px;}</style>\n";
        header += "<style>h3 {color:#000000; font-size:0.9em; font-weight:italic; padding:1px; margin:1px;}</style>\n";
        header += "<style>p  {color:" + color + "; font-size:0.9em; padding:1px; margin:0px;}</style>\n";
        header += "<style>pre  {font-size:0.8em; padding:1px; margin:0px;}</style>\n";
        header += "</head>\n";
        header += "<body>\n";
        footer += "</body></html>\n";
        setEditable(false);
        setContentType("text/html; charset=ISO-8859-1");

        append("<div style=\"text-align:center\">");
        append("h1", Constants.softname + " " + Constants.version);
        append("p", Constants.date);
        append("<p><u>" + Constants.link + "</u></p>");
        append("<p><i>" + Constants.author + "</i></p>");
        append("</div>");

        append("h2", "Reference on the method");
        append("p", ref1);

        append("h2", "Reference on OrientationJ");
        append("p", ref0);

        append("h2", "BibTeX");
        append("pre", bib0);

    }

    @Override
    public String getText() {
        Document doc = this.getDocument();
        try {
            return doc.getText(0, doc.getLength());
        }
        catch (BadLocationException e) {
            e.printStackTrace();
            return getText();
        }
    }

    public void append(String content) {
        html += content;
        setText(header + html + footer);
        setCaretPosition(0);
    }

    public void append(String tag, String content) {
        html += "<" + tag + ">" + content + "</" + tag + ">";
        setText(header + html + footer);
        setCaretPosition(0);
    }

    public JScrollPane getPane() {
        JScrollPane scroll = new JScrollPane(this, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        return scroll;
    }

    public void show(int w, int h) {
        JScrollPane scroll = new JScrollPane(this, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        scroll.setPreferredSize(new Dimension(w, h));
        JFrame frame = new JFrame();
        frame.getContentPane().add(scroll);
        frame.pack();
        frame.setVisible(true);
    }
}
