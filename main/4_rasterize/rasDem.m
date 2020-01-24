function [] = rasDem(objL,strWinPath,hW,cWin)
% Make raster grid DEM from georeferenced triangulated points

% Update waitbar
try
waitbar(str2double(cWin{1})/str2double(cWin{2}),hW, ...
    ['window ' cWin{1} ' of ' cWin{2} ': rasterizing the DEM...'])
catch
end

% Read triangulated points
mPts = objL.TriangulatedPointsGeoref;

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

% Boundaries
mB = [min(mPts(1:2,:),[],2) max(mPts(1:2,:),[],2)]';
mB(2,1:2) = mB(1,1:2) + round(diff(mB)/dRes)*dRes;
vX = mB(1):dRes:mB(2);
vY = fliplr(mB(3):dRes:mB(4));

% Interpolation parameters
sParams.blockSize = 500;
sParams.null = -9999;
sParams.connectedPixels = 0;
sParams.radius = 0;
sParams.matSource = objL.Properties.Source;
sParams.matField = 'HexagonDem';

% Interpolate to make raster DEM
objL.(sParams.matField) = [];
points2grid(mPts,vX,vY,'interp',sParams);

% Save spatial referencing structure
sR = georasterref;
sR.Lonlim = [vX(1) vX(end)] + [-dRes dRes]/2;
sR.Latlim = [vY(end) vY(1)] + [-dRes dRes]/2;
sR.RasterSize = size(objL,sParams.matField);
sR.ColumnsStartFrom = 'north';
sR.RowsStartFrom = 'west';
objL.HexagonDemSpatialRef = sR;

% Write a geotiff file
strSaveFile = strcat([strWinPath 'dems\dem_r' ...
    num2str(objL.RegionID) 'w' ...
    num2str(objL.WindowID) '.tif']);
geotiffwrite(strSaveFile,objL.HexagonDem,objL.HexagonDemSpatialRef, ...
    'CoordRefSysCode','EPSG:4326');