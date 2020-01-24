# HEXIMAP

HEXagon IMagery Automated Pipeline. MATLAB code for extracting digital elevation models (DEMs) and orthoimages from declassified HEXAGON satellite imagery (Maurer and Rupper, 2015).

## Requirements

* MATLAB version 2018a or newer

* MATLAB image processing, mapping, statistics, and optimization toolboxes.  Type `ver` in the MATLAB command window to see if you have these required toolboxes.

* [OpenCV](https://opencv.org) library

* [mexopencv](https://kyamagu.github.io/mexopencv/)

* Check [EarthExplorer](http://earthexplorer.usgs.gov/) for availability of Hexagon imagery (with minimal cloud cover) over your region of interest. Note: HEXIMAP uses an automated pipeline that involves estimating the fundamental matrix for two-view geometry. If an observed scene is nearly planar, the fundamental matrix can only be determined up to three degrees of freedom and degenerate cases can occur. So if your region of interest contains only low-relief terrain, HEXIMAP will likely fail.

## Getting started

HEXIMAP uses OpenCV APIs for things like feature detection, disparity map computation, and stereo image rectification. The software package [mexopencv](https://kyamagu.github.io/mexopencv/) acts as the interface between MATLAB and OpenCV. Detailed instructions for setting up the environment are given in the mexopencv [readme](https://github.com/kyamagu/mexopencv) and [wiki](https://github.com/kyamagu/mexopencv/wiki). 

After downloading the HEXIMAP repository, add it to your MATLAB path including all subdirectories using `addpath(genpath('/path/to/heximap'))`. Next add mexopencv to your MATLAB path using `addpath('/path/to/mexopencv')`.

## Usage

## References

* Maurer, J., & Rupper, S. (2015). Tapping into the Hexagon spy imagery database: A new automated pipeline for geomorphic change detection. ISPRS Journal of Photogrammetry and Remote Sensing, 108, 113-127.
