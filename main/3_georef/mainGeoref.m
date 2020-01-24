% Main script for georeferencing extracted Hexagon DEMs

% Waitbar handle
hW = waitbar(0/8,'please wait...');

% Check input data for errors
geoCheckInput(strRef,strShpPath,hW);

% Get mat files containing data for image windows
[cL,cR] = geoGetMatFiles(strWinPath,hW);

% User chooses window for georeferencing
iWinIdx = '';
while isempty(iWinIdx)
    iWinIdx = geoChooseWindow(cL,cR,hW);
end

% Define control points
geoControlPoints(cL{iWinIdx},strRef,strCP,hW);

% Compute transformations for initial georeferencing (using chosen window)
geoInitTrans(cL{iWinIdx},strRef,strShpPath,lVis,hW);

% Apply georeferencing transformations (from chosen window to all windows)
geoApplyTrans(cL,cR,iWinIdx,hW);

% Optimize orientation of all windows
geoOptimize(cL,iWinIdx,strRef,strShpPath,lVis,lPoly,hW);

% Update waitbar
try
waitbar(8/8,hW,'finished.')
catch
end
