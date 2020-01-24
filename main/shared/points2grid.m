function mGrid = points2grid(mPts,vX,vY,strType,varargin)
% Convert irregularly-spaced points to a grid, using sparse matrices or
% linear interpolation. For linear interpolation, can use block-processing
% to save memory.

% Make sure boundaries are correct
if vX(1) > vX(end)
    error('Horizontal coordinates must be in ascending order.')
end
if vY(1) < vY(end)
    error('Vertical coordinates must be in descending order.')
end

% Warn if user specifies sparse method and includes interp parameters
if strcmp(strType,'sparse') && ~isempty(varargin)
    warning('Input parameters will be ignored when using sparse method.')
end

% Compute resolution
dX = mode(abs(diff(vX)));
dY = mode(abs(diff(vY)));

% Make grid
switch strType    
    case 'interp'      
        mGrid = interpGrid(mPts,vX,vY,dX,dY,varargin);
    case 'sparse'
        mGrid = sparseGrid(mPts,vX,vY,dX,dY);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mGrid = interpGrid(mPts,vX,vY,dX,dY,cVarargin)

% Parameters
if ~isempty(cVarargin)
    sParams = cVarargin{1};
else
    sParams = struct();
end

% Default parameters
if ~isfield(sParams,'blockSize')
    sParams.blockSize = 500;
end
if ~isfield(sParams,'null')
    sParams.null = NaN;
end

% Check if user wants to mask holes
lMask = false;
if isfield(sParams,'connectedPixels') && isfield(sParams,'radius')
if isnumeric(sParams.connectedPixels) && isnumeric(sParams.radius)
    lMask = true;
end
end

% Get matfile path and field if user specifies
lMat = false;
if isfield(sParams,'matSource') && isfield(sParams,'matField')
if ischar(sParams.matSource) && ischar(sParams.matField)
    sMatInfo.file = sParams.matSource;
    sMatInfo.field = sParams.matField;
    try
        sMatInfo.startIndex = sParams.matStartIndex;
    catch
        sMatInfo.startIndex = [1 1];
    end
    lMat = true;
end
end

% Edge buffer size
dB = mean([dX dY])*10;

% Size of spatial referencing vectors
vSz = [length(vY) length(vX)];

% Function handle for processing blocks
hFun = @(x) makeBlock(x,mPts,sParams,vX,vY,dX,dY,dB,lMask);

% Call block processing function
if lMat
    mGrid = blockProcess(vSz,sParams.blockSize,'double',hFun,sMatInfo);
else
    mGrid = blockProcess(vSz,sParams.blockSize,'double',hFun);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mBlock = makeBlock(sInput,mPts,sParams,vX,vY,dX,dY,dB,lMask)

% Get input from block processing function
vIdxX = sInput.indexX;
vIdxY = sInput.indexY;
iX = sInput.countX;
iY = sInput.countY;

% Define current block
vBlkX = [vX(vIdxX(iX)) vX(vIdxX(iX+1))];
vBlkY = [vY(vIdxY(iY)) vY(vIdxY(iY+1))];
    
% Define inliers
lIn = (mPts(1,:) >= vBlkX(1)-dB & mPts(1,:) <= vBlkX(2)+dB) & ...
      (mPts(2,:) >= vBlkY(2)-dB & mPts(2,:) <= vBlkY(1)+dB);

% Skip current block if empty
if sum(lIn) == 0
    mBlock = [];
    return
end

% Compute interpolant
F = scatteredInterpolant(mPts(1,lIn)',mPts(2,lIn)',mPts(3,lIn)', ...
    'linear','none');

% Spatial referencing for current block
vIdxXb = vIdxX(iX):vIdxX(iX+1);
vIdxYb = vIdxY(iY):vIdxY(iY+1);
[mXb,mYb] = meshgrid(vX(vIdxXb),vY(vIdxYb));

% Sample interpolant
mBlock = F(mXb,mYb);

% Mask holes
if lMask
    lM = isnan(sparseGrid(mPts(:,lIn), ...
         sort(vBlkX,'ascend'),sort(vBlkY,'descend'),dX,dY));
    lM = bwareaopen(lM,sParams.connectedPixels);
    lM = imdilate(lM,strel('disk',sParams.radius));
    mBlock(lM) = NaN;
end

% Replace NaNs with user-specified null value
mBlock(isnan(mBlock)) = sParams.null;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mGrid = sparseGrid(mPts,vX,vY,dX,dY)

% Flip y-axis direction for matrix coordinate system
vY = flipud(vY(:))';
    
% Assign matrix indices to each point
mPts(1,:) = round((mPts(1,:) - vX(1))/dX + 1);
mPts(2,:) = round((mPts(2,:) - vY(1))/dY + 1);

% Matrix size
vSz(1) = round((vX(end) - vX(1))/dX + 1);
vSz(2) = round((vY(end) - vY(1))/dY + 1);

% Remove any points outside boundaries
lM = mPts(2,:) > vSz(2) | mPts(2,:) < 1 | ...
     mPts(1,:) > vSz(1) | mPts(1,:) < 1;
mPts(:,lM) = [];

% Make grid
mGrid = sparse(mPts(2,:),mPts(1,:),mPts(3,:),vSz(2),vSz(1));
mCount = sparse(mPts(2,:),mPts(1,:),ones(1,length(mPts)),vSz(2),vSz(1));
mGrid = full(mGrid ./ mCount);

% Flip y-axis back to original orientation
mGrid = flipud(mGrid);

end
end
