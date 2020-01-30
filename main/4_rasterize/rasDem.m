function [] = rasDem(objL,strWinPath,lClean,lMed,lDen,iGap,iSpec,iMed, ...
    iDenT,iDenN,hW,cWin)
% Make raster grid DEM from georeferenced triangulated points

% Update waitbar
try
waitbar(str2double(cWin{1})/str2double(cWin{2}),hW, ...
    ['window ' cWin{1} ' of ' cWin{2} ': rasterizing the DEM...'])
catch
end

% Read triangulated points
mPts = objL.TriangulatedPointsGeoref;
mPts(:,any(isnan(mPts),1)) = [];

% Compute reasonable resolution in geographic coordinate system, given the
% spacing of the triangulated points
sGeo = objL.GeorefInfo;
sT = sGeo.Initial.Triangulated2WorldTransform;
vRes = abs(diff(round(ll2utm(mPts',sT.zone,sT.hemi))))';
[vCount,vRes] = hist(vRes(1,:),unique(vRes(1,:)));
[~,vIdx] = sort(vCount,'descend');
dResM = max([mean(vRes(vIdx(1:3))) 1]);
clear vRes vCount vIdx
vC = ll2utm(mean(mPts(1:2,:),2)',sT.zone,sT.hemi);
dRes = mean(abs(diff(utm2ll([vC;vC + dResM],sT.zone,sT.hemi))));
dRes = dRes*2;
dResM = dResM*2;

% Boundaries
mB = [min(mPts(1:2,:),[],2) max(mPts(1:2,:),[],2)]';
mB(2,1:2) = mB(1,1:2) + round(diff(mB)/dRes)*dRes;
vX = mB(1):dRes:mB(2);
vY = fliplr(mB(3):dRes:mB(4));

% Interpolation parameters
sParams.blockSize = 500;
sParams.null = NaN;
sParams.connectedPixels = 0;
sParams.radius = 0;

% Interpolate to make raster DEM
mDem = points2grid(mPts,vX,vY,'interp',sParams);

% Clean up DEM
mDem = rasClean(mDem,dResM,lClean,iGap,iSpec);

% Median filter and mesh denoise 
mDem = rasSmooth(vX,vY,mDem,lMed,lDen,iMed,iDenT,iDenN);

% Set nodata value and numeric class
mDem(isnan(mDem)) = -32768;
mDem = int16(mDem);

% Spatial referencing structure
sR = georasterref;
sR.Lonlim = [vX(1) vX(end)] + [-dRes dRes]/2;
sR.Latlim = [vY(end) vY(1)] + [-dRes dRes]/2;
sR.RasterSize = size(mDem);
sR.ColumnsStartFrom = 'north';
sR.RowsStartFrom = 'west';

% Write a geotiff file
strSaveFile = strcat([strWinPath 'dems\dem_r' ...
    num2str(objL.RegionID) '_w' ...
    num2str(objL.WindowID) '.tif']);
geotiffwrite(strSaveFile,mDem,sR,'CoordRefSysCode','EPSG:4326');

% Save output
objL.HexagonDem = mDem;
objL.HexagonDemSpatialRef = sR;