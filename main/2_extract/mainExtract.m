% Main script for extracting DEMs from Hexagon images

% Waitbar handle
hW = waitbar(0/8,'please wait...');

% Define Hexagon mat files
cM = cellfun(@(x) matfile(strcat([strHexPath x]),'Writable',true), ...
    cHexFile, 'UniformOutput',0);

% Make sure user selected exactly two Hexagon mat files
if numel(cM) ~= 2
    error('You must select exactly two Hexagon .mat files (stereo pair).')
end

% User provides image corner points (from earthexplorer.usgs.gov) for
% roughly georeferencing and aligning the images
cellfun(@(x,y) extControlPoints(x,y,hW),cM,cHexFile,'UniformOutput',0);

% Sort the Hexagon images to establish left and right images when computing
% stereo disparity maps
[cHexFile,cM] = extSortImages(cHexFile,cM,hW);

% Compute fundamental matrix, relative pose matrices, and rough
% homographies using full images. More accurate matrices are computed
% later for each processing window.
extInitTrans(cM,strHexPath,hW);

% User selects regions of interest
cWindow = extChooseWindows(cM,hW);

% Rectify the stereo images, compute disparity maps
extDisparityLoop(cM,strHexPath,cWindow,strRes,iBlkSz,hW);

% Refine camera orientations using bundle adjustment
extBundleAdjustLoop(strHexPath,hW);

% Triangulate the points
extTriangulateLoop(strHexPath,hW);

% Update waitbar
try
waitbar(8/8,hW,'finished.')
catch
end
