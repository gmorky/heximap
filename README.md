# HEXIMAP

HEXagon IMagery Automated Pipeline. MATLAB code for extracting digital elevation models (DEMs) and orthoimages from declassified KH-9 Hexagon satellite imagery (Maurer and Rupper, 2015; Maurer et al., 2019).

## Requirements

* MATLAB version 2018a or newer

* MATLAB image processing, mapping, statistics, and optimization toolboxes. Enter `ver` in the MATLAB command window to see if you have them installed.

* [OpenCV](https://opencv.org) library

* [mexopencv](https://kyamagu.github.io/mexopencv/)

* Digital scans of declassified Hexagon imagery provided by the United States Geological Survey (USGS). Check [EarthExplorer](http://earthexplorer.usgs.gov/)  (Data Sets >> Declassified Data >> Declass 2) for availability of cloud-free Hexagon images over your region of interest. The USGS can scan at 7 μm (higher resolution) or 14 μm (lower resolution). Currently HEXIMAP only supports 7 μm scans of film from the KH-9 Hexagon "Lower Resolution Mapping Camera".

## Tips

* Any external data used as input must be georeferenced in the WGS84 geographic coordinate system, with elevations specified in meters.

* HEXIMAP uses computer vision algorithms to estimate the fundamental matrix for two-view geometry. If an observed scene is nearly planar, the fundamental matrix can only be determined up to three degrees of freedom and degenerate cases can occur. If your region of interest only contains very flat terrain, HEXIMAP will likely fail.

* The current automatic georeferencing code can struggle to achieve accurate results in some situations. When possible, users should check the accuracy of HEXIMAP output using independently georeferenced data.

* Oftentimes Hexagon film has poor image contrast over bright or featureless terrain; particularly for snow, shadows, and water. The stereo matching algorithm fails in these situations, and output DEMs will often have data gaps or inaccurate elevations. Some manual cleanup and removal of erroneous DEM pixels from HEXIMAP output is usually required, and caution should be used when analyzing and interpreting the data over such regions.

## Installation

HEXIMAP relies on OpenCV APIs for things like feature detection, disparity map computation, and stereo image rectification. The software package [mexopencv](https://kyamagu.github.io/mexopencv/) acts as the interface between MATLAB and OpenCV. Detailed instructions for setting up the environment are given in the mexopencv [readme](https://github.com/kyamagu/mexopencv) and [wiki](https://github.com/kyamagu/mexopencv/wiki). 

After downloading the HEXIMAP repository, add it to your MATLAB path including all subdirectories as `addpath(genpath('/path/to/heximap'))`. Also add mexopencv as `addpath('/path/to/mexopencv')`.

## Usage

*Any external data used as input must be georeferenced in the WGS84 geographic coordinate system, with elevations specified in meters.* There are four stages in the HEXIMAP workflow. Each stage has a corresponding GUI where inputs and parameters can be specified before running:

#### Stitch 

Enter `stitch` in the MATLAB command window. Hexagon image halves are stitched together using feature extraction and matching. This is required because film strips are digitally scanned by the USGS as left and right halves with some overlap.

Input: 

* Select at least two Hexagon geotiff files (left half and right half) downloaded from [EarthExplorer](http://earthexplorer.usgs.gov/). Left image halves have "a" in the filename, while right image halves have "b" in the filename. For example, a Hexagon image with ID "DZB1210-500177L001001" would include "DZB1210-500177L001001_a.tif" and "DZB1210-500177L001001_b.tif". *Sometimes these labels are reversed when downloaded from the EarthExplorer website, so make sure to rename the image files if incorrect (by visual inspection of the black border around the image edges) before proceeding.* You can select more than two images as input (i.e. multiple “a-b” pairs) if you wish. Just make sure corresponding pairs have correct filenames as described.

Additional steps:

* Before stitching the image halves, HEXIMAP requires approximate image corner locations. It will attempt to find them automatically, but may fail to do so. For each corner a figure will open showing the detected corner, and a dialogue will ask to confirm whether it was detected correctly. If the small orange marker appears to be in the correct location, click “yes”. Otherwise, click “no” and manually specify the corner by clicking on it when prompted.  

* After all corners have been located, HEXIMAP will stitch the image halves together and save the full Hexagon images as \*.mat files in the same folder as the geotiff files from EarthExplorer.

#### Extract

Enter `extract` in the MATLAB command window. DEMs are extracted from a pair of Hexagon images over a region of interest. This is accomplished by estimating the fundamental matrix and relative camera pose matrices between an image pair via feature matching, then applying homography transforms for stereo rectification. Disparity maps are computed using semi-global block matching (SGBM) (Hirschmüller, 2008), relative camera orientations are refined using a bundle adjustment procedure, and DEM points are triangulated in the camera coordinate system.

Input: 

* Select exactly two \*.mat files containing the stitched Hexagon images (output from the "stitch" stage). The two files must be a pair with stereo overlap, acquired as the satellite proceeded along its orbital trajectory over a region. Usually the USGS filenames correspond to stereo pairs or triplets. For example, if three Hexagon images downloaded over a region of interest were named “DZB1216-500488L004001” DZB1216-500488L005001”, and “DZB1216-500488L006001”, the numbers 4, 5, and 6 in each respective filename would correspond to the order the images were acquired by the camera. You would want to select either 4 and 5, or 5 and 6, but not 4 and 6. This is because the baseline between 4 and 6 would likely be too wide for stereo matching, especially for mountainous terrain.

* Choose the desired image resolution for stereo matching. In many situations, the “1/2” setting provides satisfactory results. The "Full" setting is much more computationally expensive and often results in DEMs with more speckles and data gaps due to noise in the Hexagon imagery and the nature of the SGBM algorithm.

* SGBM block size: This is the size of the block (in pixels) used for matching between stereo images when creating the disparity map. It should be an odd number, usually ranging between 3 to 11.

Additional steps:

* HEXIMAP requires approximate ground control points for the image corners for initial georeferencing. After clicking "run", a dialogue will open where you can enter ground control points for the image corners. These can be found on [EarthExplorer](http://earthexplorer.usgs.gov/), by clicking on the “show metadata and browse” icon for an image of interest. Be sure to enter the correct latitude and longitude values in decimal degrees for each image corner. Note: you must do this twice - once for each image in the pair. The current image is given by the label at the top of the dialogue box.

* Next a figure will appear showing a lower-resolution version of one of the Hexagon images, and a draggable and resizable box (also referred to as “window” in the georeferencing section below) allowing you to select image regions to process. Refer to Google Earth or similar program if an external visual reference is needed (the Hexagon image will be sideways, with north-south direction running horizontally). To select a region, double click inside the window box. A dialogue will then appear asking if you wish to specify another region of interest. Select “yes” to continue selecting multiple regions. When you select "no", the DEM extraction process will begin. *Important: Each time you double click inside the box to select a region, that region will be georeferenced independently. You may need to experiment to find an optimal region size for accurate georeferencing. Selecting a single region any larger than 2 boxes across and 2 boxes wide is usually not recommended.* Also note the vertical orange dotted line, which indicates the boundary for the region of stereo overlap between the two Hexagon images. Only this overlap region is viable for processing in HEXIMAP. If needed, you can use the zoom and pan tools located at the top of the figure to explore the image.

#### Georeference

Enter `georef` in the MATLAB command window. Extracted DEMs are georeferenced in a standard world coordinate system (WGS84). HEXIMAP attempts to shift, rotate, and scale the Hexagon DEMs to match an existing georeferenced DEM (SRTM usually works well) using terrain feature matching, an iterative closest point (ICP) algorithm, and nonlinear optimization.

Input:

* Select a reference DEM. This must be a single geotiff file in standard WGS84 geographic coordinate system.  For most regions, an open source global DEM such as SRTM is ideal.

* Select a folder containing ESRI shapefiles (\*.shp) of polygons enclosing any unstable terrain in your region of interest such as glaciers, landslides, coastal erosion, or any other areas where significant geomorphic change has occurred.  These polygons do not need to be precise, but should mask approximate areas of unstable terrain to avoid elevation biases during georeferencing. The input shapefile(s) must use the WGS84 geographic coordinate system. If you do not suspect any unstable terrain exists in your region of interest, you can leave this input blank.

* Select the top folder containing all extracted DEM files from the previous stage.

* Choose a control point selection method. The “automatic” setting will attempt to find terrain feature matches between the Hexagon DEMs and the reference DEM. Sometimes this automated process fails, and you will need to rerun the code with the “manual” method selected. Then ground control points will need to be input manually.

* Elevation difference outlier threshold: Threshold for which (absolute) elevation differences between the Hexagon DEM and reference DEM are considered outliers. Any Hexagon DEM values with a difference greater than this threshold are removed. This can be useful for filtering out erroneous elevations caused by clouds or poor image contrast.

* Show animation during optimization: This opens a figure to illustrate progress of georeferencing during the optimization.

* Use polynomial surface correction: This applies a 3rd order polynomial surface correction to the Hexagon DEMs, which can sometimes help reduce georeferencing errors due to film distortions in the Hexagon imagery. However, this option can sometimes introduce new errors into the DEMs, particularly in regions where unstable terrain has been masked along DEM edges. Use with caution.

* Perform additional optimization of individual windows: It is possible to for a user to define a region spanning multiple “windows” (also referred to as “boxes” in the DEM extraction section above). If this option is selected, each window within a single region will receive an additional optimization run after the primary optimization of the larger “region”. This can sometimes improve the overall georeferencing accuracy. Note: only windows with greater than 33% terrain data coverage will receive the additional optimization step.

Additional steps:

* After clicking “run”, a figure will appear showing a lower-resolution version of one of the Hexagon images, along with all the labelled windows (orange boxes) for which DEMs have been extracted. For regions containing multiple windows, you’ll notice that the boxes overlap one another to ensure minimal data gaps along edges.

* A dialogue will also appear, where you can enter the number ID of the window (box) to use for initial georeferencing. This should be a window with sufficient image contrast, minimal cloud cover, some vertical relief in the terrain, a good disparity map, and should not contain a significant amount of unstable terrain. Enter the number ID of a window you suspect is of high quality, then another figure will open containing a closer view of the chosen window with a disparity map overlay. Cool colors are further from the camera (lower elevations) while warm colors are closer to the camera (higher elevations). If the disparity map appears to be of good quality without too many large gaps, select “Use this region” at the bottom of the figure. Otherwise, select “Choose new region”. You can also click on the disparity map to toggle the visibility. After selecting a region, HEXIMAP will automatically georeference the remaining regions.

#### Rasterize

Enter `rasterize` in the MATLAB command window. Georeferenced DEMs and orthoimages are rasterized and exported as geotiff files. Before exporting, the raster DEMs can (optionally) be cleaned up, spatially filtered, and denoised as described below. Exported rasters and orthoimages are saved in folders “/dems”, and “/images”, respectively. Any data gaps in the DEMs are assigned a NoData value of -32768.

Input:

* Select top folder containing the extracted Hexagon DEM files (should be the same folder from the previous stages).

* Cleanup DEMs: Additional cleanup steps on the DEMs before exporting. This includes interpolating small gaps, a slight dilation of larger gaps to remove erroneous pixels along gap edges, and removal of isolated “speckles” i.e. small islands of pixels which are completely surrounded by NoData values (these speckles are often erroneous elevations).

	* Gap threshold: Any data gaps smaller than this value (sq. meters) will be filled using bilinear interpolation.

  * Speckle threshold: Any isolated speckles smaller than this value (sq. meters) will be removed.

* Filter DEMs:  2-D median filtering of the DEMs.

  * Window size: Size of the sliding window for the median filter.

* Denoise DEMs: This final processing step can help reduce noise in the DEMs (caused by noise in the Hexagon imagery and issues inherent to the stereo matching algorithm) while still preserving terrain features (Sun et al., 2007).

  * Param t: The threshold parameter in Sun et al. (2007) which determines how many surface normals will be used in the locally weighted averaging operation. Must be between 0 and 1. A larger value for t better preserves sharp terrain features (such as ridge crests) but allows more noise to pass through.

  * Param n: This corresponds to the number of surface normal updating operations. A higher value leads to greater smoothing of the DEMs.

## References

* Maurer, J., & Rupper, S. (2015). Tapping into the Hexagon spy imagery database: A new automated pipeline for geomorphic change detection. ISPRS Journal of Photogrammetry and Remote Sensing, 108, 113-127.

* Maurer, J. M., Schaefer, J. M., Rupper, S., & Corley, A. (2019). Acceleration of ice loss across the Himalayas over the past 40 years. Science advances, 5(6), eaav7266.

* Hirschmüller, H. (2008). Stereo processing by semiglobal matching and mutual information. IEEE transactions on pattern analysis and machine intelligence, 30(2), 328-341.

* Sun, X., Rosin, P. L., Martin, R., & Langbein, F. (2007). Fast and effective feature-preserving mesh denoising. IEEE transactions on visualization and computer graphics, 13(5), 925-938.
