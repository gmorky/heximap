function [mData,vX,vY] = readGeotiffRegion(mWin,strPath,varargin)
% Read pixel region of a geotiff file

% Check for valid pixel region
if any(size(mWin) ~= 2)
    error('Invalid pixel region.')
end

% Apply edge buffer
mWin = sort(mWin,1);
if ~isempty(varargin)
    dBuffPct = varargin{1};
    vB = diff(mWin,1)*dBuffPct/100;
    mWin = mWin + [-vB;vB];
end

% Get spatial referencing info
sInfo = geotiffinfo(strPath);
sR = sInfo.SpatialRef;

% Make spatial referencing vectors
[vX,vY] = makeSpatialRefVecs(sR,'full');

% Find the closest data pixels to the edges of the specified window
vIdxX = sort(knnsearch(vX',mWin(:,1)));
vIdxY = sort(knnsearch(vY',mWin(:,2)));
vX = vX(vIdxX(1):vIdxX(2));
vY = vY(vIdxY(1):vIdxY(2));

% Error check
if diff(vIdxX)<1 || diff(vIdxY)<1
   error(['Empty dataset. Check coordinate systems of the ' ...
   'input data, and make sure the specified window and geotiff file ' ...
   'have spatial overlap.'])
end

% Read the data
mData = imread(strPath,'PixelRegion', ...
    {[vIdxY(1) vIdxY(2)] [vIdxX(1) vIdxX(2)]});
