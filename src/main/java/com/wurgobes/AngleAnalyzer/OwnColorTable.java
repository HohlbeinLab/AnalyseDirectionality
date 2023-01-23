package com.wurgobes.AngleAnalyzer;

import net.imagej.lut.LUTService;
import net.imglib2.display.ColorTable;
import org.scijava.command.Command;
import org.scijava.plugin.Plugin;

import java.awt.*;

/**
 * From <a href="https://forum.image.sc/t/plugins-java-how-to-get-a-color-from-current-specific-lut-using-an-index/4685/10">...</a>
 * @author alex.vergara
 */
@Plugin(type = Command.class)
public class OwnColorTable {
    private Color[] CT;

    private final LUTService ls;
    private boolean shownFailure = false;

    private String currentLUT;

    public OwnColorTable(LUTService ls){
        this.ls = ls;
    }

    public void setLut(String colormap) {
        if(ls == null){
            backupLUT();
        } else {
            try {
                ColorTable ct = ls.loadLUT(ls.findLUTs().get(colormap));
                int size = ct.getLength();
                CT = new Color[size];
                for (int i = 0; i < size; i++) {
                    CT[i] = new Color(ct.get(ColorTable.RED, i), ct.get(ColorTable.GREEN, i), ct.get(ColorTable.BLUE, i));
                }
                currentLUT = colormap;
            }
            catch (Exception e){
                //e.printStackTrace();
                backupLUT();
            }
        }

    }

    private void backupLUT(){
        currentLUT = "fire.lut";
        if(!shownFailure) {
            System.out.println("LUT Service failed, using backup LUT: fire.lut");
            shownFailure = true;
        }
        int[] r = {0,0,1,25,49,73,98,122,146,162,173,184,195,207,217,229,240,252,255,255,255,255,255,255,255,255,255,255,255,255,255,255};
        int[] g = {0,0,0,0,0,0,0,0,0,0,0,0,0,14,35,57,79,101,117,133,147,161,175,190,205,219,234,248,255,255,255,255};
        int[] b = {0,61,96,130,165,192,220,227,210,181,151,122,93,64,35,5,0,0,0,0,0,0,0,0,0,0,0,35,98,160,223,255};
        CT = new Color[r.length];
        for (int i=0; i<r.length; i++) {
            CT[i] = new Color(r[i], g[i], b[i]);
        }
    }


    public Color getColor(float value, double min, double max){
        int v = (int)(Math.max(Math.min((value - min)/(max - min), 1), 0) * (CT.length-1)); //remap from 0 to number of lut colors

        return CT[v];
    }

    public String[] getLuts(){
        return ls.findLUTs().keySet().toArray(new String[0]);
    }

    public String getCurrentLUT(){return currentLUT;}

    public Color getColor(float zValue, double[] lutRange) {
        return getColor(zValue, lutRange[0], lutRange[1]);
    }
}