function [] = rasOrtho(objL,strWinPath,hW,cWin)
% Make raster grid orthoimage

% Update waitbar
try
waitbar(str2double(cWin{1})/str2double(cWin{2}),hW, ...
    ['window ' cWin{1} ' of ' cWin{2} ': rasterizing the orthoimage...'])
catch
end

% Get georeferencing info
sGeo = objL.GeorefInfo;
sT = sGeo.Initial.Triangulated2WorldTransform;
iZ = sT.zone;
strH = sT.hemi;

% Make spatial referencing vectors for Hexagon DEM
sR = objL.HexagonDemSpatialRef;
[vLonH,vLatH] = makeSpatialRefVecs(sR,'full');

% Make spatial referencing vectors for new DEM
dRes = 10;
vC = ll2utm([median(vLonH) median(vLatH)],iZ,strH);
vC = utm2ll([vC;vC+dRes],iZ,strH);
dRes = mean(abs(diff(vC)));
vLon = vLonH(1):dRes:vLonH(end);
vLat = vLatH(1):-dRes:vLatH(end);

% Make new DEM
[mLon,mLat] = meshgrid(vLon,vLat);
[mLonH,mLatH] = meshgrid(vLonH,vLatH);
mDem = objL.HexagonDem; mDem(mDem == -9999) = NaN;
mDem = interp2(mLonH,mLatH,mDem,mLon,mLat);
clear mLonH mLatH

% Fill holes with reference DEM
[mLonR,mLatR] = meshgrid(objL.ReferenceLon,objL.ReferenceLat);
lHoles = isnan(mDem);
mDem(lHoles) = interp2( ...
    mLonR,mLatR,objL.ReferenceDem,mLon(lHoles),mLat(lHoles));
clear mLonR mLatR lHoles

% Interpolate any remaining holes (where the reference DEM also had no
% data)
mDem = inpaint_nans(mDem);

% Make DEM points
vSz = size(mDem);
mPts = [mLon(:) mLat(:) mDem(:) ones(numel(mLon),1)]';
clear mLon mLat mDem

% Transform points from WGS84 latitude and longitude to UTM
mPts = ll2utm(mPts',iZ,strH)';

% Undo final optimizations and alignments
cOutput = sGeo.Final.OptimizationOutput;
cT = sGeo.Final.AlignmentOutput;
for iT = numel(cOutput):-1:1
    sOutput = cOutput{iT};
    sOutput.direction = 'inverse';
    mPts = transformUsingSolverVar(mPts,sOutput);
    mPts = [cT{iT};0 0 0 1] \ mPts;
end

% Undo curvature correction
mPts = curvatureCorrection(mPts,iZ,strH,'inverse',sT.curvature);

% Undo initial optimization and alignment
sOutput = sGeo.Initial.OptimizationOutput;
sOutput.direction = 'inverse';
mPts = transformUsingSolverVar(mPts,sOutput);
mPts = [sGeo.Initial.AlignmentOutput;0 0 0 1] \ mPts;

% Undo initial transformations. Now the reference DEM points are in the
% original coordinate system of the triangulated points.
mPts = sT.trans \ mPts;
if ~strcmp(sGeo.Initial.WindowTransform,'none')
    mPts = sGeo.Initial.WindowTransform \ mPts;
end

% Project onto full image plane
mPts = objL.IntrinsicMatrix * objL.PoseMatrix * mPts;
mPts(1,:) = mPts(1,:) ./ mPts(3,:);
mPts(2,:) = mPts(2,:) ./ mPts(3,:);
mPts(3,:) = [];

% Put into image window coordinates
mHexWin = objL.Window;
mPts(1,:) = mPts(1,:) - mHexWin(1) + 1;
mPts(2,:) = mPts(2,:) - mHexWin(3) + 1;

% Interpolate to make raster orthoimage
objL.HexagonImage = ...
    uint8(reshape(interp2(double(objL.Image),mPts(1,:),mPts(2,:)),vSz));

% Save spatial referencing structure
sR = georasterref;
sR.Lonlim = [vLonH(1) vLonH(end)] + [-dRes dRes]/2;
sR.Latlim = [vLatH(end) vLatH(1)] + [-dRes dRes]/2;
sR.RasterSize = vSz;
sR.ColumnsStartFrom = 'north';
sR.RowsStartFrom = 'west';
objL.HexagonImageSpatialRef = sR;

%Write a geotiff file
strSaveFile = strcat([strWinPath 'images\image_r' ...
    num2str(objL.RegionID) 'w' ...
    num2str(objL.WindowID) '.tif']);
geotiffwrite(strSaveFile,objL.HexagonImage,objL.HexagonImageSpatialRef, ...
    'CoordRefSysCode','EPSG:4326');
