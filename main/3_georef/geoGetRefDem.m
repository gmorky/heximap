function [mDem,vLon,vLat] = geoGetRefDem(mWin,strRef)

% Read the geotiff file
try
    [mDem,vLon,vLat] = readGeotiffRegion(mWin,strRef,10);
catch objExc
    warning(objExc.message)
    error('Failed to read reference DEM.')
end

% Convert to double and exclude any null values
mDem = double(mDem);
mDem(mDem > 9000 | mDem < -500) = NaN;
