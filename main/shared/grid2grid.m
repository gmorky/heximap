function mGrid = grid2grid(sourceGridInfo,vX,vY,varargin)
% Resample a regularly-spaced grid

% Check input parameters
if isempty(varargin)
    sParams = struct();
else
    sParams = varargin{1};
end

% Make sure boundaries are correct
if vX(1) > vX(end)
    error('Horizontal coordinates must be in ascending order.')
end
if vY(1) < vY(end)
    error('Vertical coordinates must be in descending order.')
end
 
% Default block size
if ~isfield(sParams,'blockSize')
    sParams.blockSize = 500;    
end

% Default null value
if ~isfield(sParams,'nullVal')
    sParams.nullVal = NaN;
end

% Default tiff layer
if ~isfield(sParams,'tiffLayer')
    sParams.tiffLayer = 1;
end

% Check transformation
if ~isfield(sParams,'transform')
    sParams.transform = [];
else
    if isstruct(sParams.transform)
        if all(~strcmp(sParams.transform.direction, ...
                {'ll2utm','utm2ll','ll2ps','ps2ll','utm2ps','ps2utm'})) 
            error('Invalid transformation direction.');
        end
    elseif isnumeric(sParams.transform)
        if any(size(sParams.transform) ~= [4 4])
            error('Invalid transformation matrix.')
        end
    else
    end
end

% Default interpolation method
if ~isfield(sParams,'interp')
    sParams.interp = 'linear';
end

% Get input grid info (for reading grid)
if ischar(sourceGridInfo)
    sInfo = geotiffinfo(sourceGridInfo);
elseif isstruct(sourceGridInfo)
    sInfo = sourceGridInfo;
else
    error('Invalid source grid input.')
end

% Get mat file info if user specifies (for saving grid)
lMat = false;
if isfield(sParams,'matSource') && isfield(sParams,'matField')
if ischar(sParams.matSource) && ischar(sParams.matField)
    sMatSaveInfo.file = sParams.matSource;
    sMatSaveInfo.field = sParams.matField;
    try
        sMatSaveInfo.startIndex = sParams.matStartIndex;
    catch
        sMatSaveInfo.startIndex = [1 1];
    end
    lMat = true;
end
end

% Edge buffer size in pixels
iB = 100;

% Size of spatial referencing vectors
vSz = [length(vY) length(vX)];

% Function handle for processing blocks
hFun = @(x) makeBlock(x,sInfo,vX,vY,iB,sParams.nullVal, ...
    sParams.tiffLayer,sParams.transform,sParams.interp);

% Call block processing function
warning('off','MATLAB:wrongBlockSize');
if lMat
    mGrid = blockProcess(vSz,sParams.blockSize,'double',hFun,sMatSaveInfo);
else
    mGrid = blockProcess(vSz,sParams.blockSize,'double',hFun);
end
warning('on','MATLAB:wrongBlockSize');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mBlock = makeBlock(sInput,sInfo,vX,vY,iB,nullVal,iL,transform, ...
        strInterp)
    
% Get input from block processing function
vIdxX = sInput.indexX;
vIdxY = sInput.indexY;
iX = sInput.countX;
iY = sInput.countY;

% Spatial referencing grids for new image block
vBlkX = vX(vIdxX(iX):vIdxX(iX+1));
vBlkY = vY(vIdxY(iY):vIdxY(iY+1));
[mX,mY] = meshgrid(vBlkX,vBlkY);

% Apply transformation
if ~isempty(transform)
    if isstruct(transform)
        iZ = transform.zone;
        strH = transform.hemi;
        switch transform.direction
            case 'll2utm'
                mPts = utm2ll([mX(:) mY(:)],iZ,strH)';
            case 'utm2ll'
                mPts = ll2utm([mX(:) mY(:)],iZ,strH)';            
            case 'll2ps'
                mPts = ps2ll([mX(:) mY(:)],transform.polar)';
            case 'ps2ll'
                mPts = ll2ps([mX(:) mY(:)],transform.polar)';             
            case 'utm2ps'
                mPts = ps2ll([mX(:) mY(:)],transform.polar);
                mPts = ll2utm(mPts,iZ,strH)';
            case 'ps2utm' 
                mPts = utm2ll([mX(:) mY(:)],iZ,strH);
                mPts = ll2ps(mPts,transform.polar)'; 
        end        
    else
        mT = transform;
        mPts = mT \ [mX(:) mY(:) ones(numel(mX),2)]';
    end
    mX = reshape(mPts(1,:),size(mX));
    mY = reshape(mPts(2,:),size(mY));
end
clear mPts

% Spatial referencing vectors for full transform image
sR = sInfo.SpatialRef;
[vXt,vYt] = makeSpatialRefVecs(sR,'full');

% Find overlapping region of transform image
vIdxX = knnsearch(vXt',[min(mX(:)) max(mX(:))]')';
vIdxY = knnsearch(vYt',[max(mY(:)) min(mY(:))]')';

% Return if no overlap
if diff(vIdxX) < 2 || diff(vIdxY) < 2
    mBlock = [];
    return
end

% Add buffer
vIdxX = vIdxX + [-iB iB];
vIdxY = vIdxY + [-iB iB];

% Make sure block doesn't extend past boundaries of transform image
vIdxX(vIdxX < 1) = 1;
vIdxY(vIdxY < 1) = 1;
vIdxX(2) = min([vIdxX(2) sInfo.Width]);
vIdxY(2) = min([vIdxY(2) sInfo.Height]);

% Spatial referencing grids for transform image
[mXt,mYt] = meshgrid(vXt(vIdxX(1):vIdxX(end)),vYt(vIdxY(1):vIdxY(end)));

% Read transform image block
try
    mBlock = imread(sInfo.Filename,'pixelregion',{vIdxY,vIdxX});
catch
    objM = matfile(sInfo.matSource);
    eval(['mBlock = objM.' sInfo.matField ...
        '(vIdxY(1):vIdxY(end),vIdxX(1):vIdxX(end));'])
end

% Find null pixels
if ischar(nullVal)
    switch nullVal
        case 'dem'
            lM = mBlock < -500 | mBlock > 9000;
        otherwise
            error('Invalid null value string.')
    end
elseif isnumeric(nullVal)
    lM = false(size(mBlock));
    for i = 1:numel(nullVal)
        lM = lM | mBlock == nullVal(i);
    end
else
    error('Invalid null value type.')
end

% Convert to double and set null pixels to NaN
mBlock = double(mBlock);
mBlock(lM) = NaN;

% Resample
mBlock = interp2(mXt,mYt,mBlock(:,:,iL),mX,mY,strInterp,NaN);

end
end

