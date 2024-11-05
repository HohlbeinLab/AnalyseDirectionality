# Angle Analyzer

This package encompasses both the ImageJ plugin RFT, and the scripts used to further process the data.

## ImageJ
To automatically use the latest version, you can subscribe to the Hohlbein update site in ImageJ.
Otherwise, download the latest version from [link](https://github.com/HohlbeinLab/AnalyseDirectionality/releases) and place it in the plugins folder in your FIJI installation folder and restart the application.

To use the plugin, open an image and click Plugins > Angle Analyzer > Analyze Angles
This opens the following window:

(insert image here)

It contains the following settings:
 1. 
 2. 

### Macro

The plugin can be run from a macro using the following syntax:

```
run("Analyze Angles", "buffer=0 window=101 overlap=0.75 path='path\\to\\folder with spaces'");
```	

All arguments follow the format of `keyword=value` using the period as the separator. When a path contains spaces, use single quote ticks to ensure the proper path is captured

The arguments are as follows:
 - buffer
 - window
 - cutoff
 - overlap
 - intensity_cutoff
 - start
 - end
 - step
 - scanning_range
 - scan_save

 - vector_overlay
 - vector_length
 - vector_width

 - path
 - order