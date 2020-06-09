function mGrid = blockProcess(vSz,dBlkSz,strClass,hFun,varargin)

% Get mat file info
lMat = false;
if nargin == 5
    sMatInfo = varargin{1};
    strFile = sMatInfo.file;
    strField = sMatInfo.field;
    try
        vStartIdx = sMatInfo.startIndex;
    catch
        vStartIdx = [1 1];
    end
    objM = matfile(strFile,'Writable',true);
    lMat = true;
elseif nargin > 5
    error('Too many input arguments.')
end

% Block vectors
iNumBlkX = ceil(vSz(2)/dBlkSz);
iNumBlkY = ceil(vSz(1)/dBlkSz);
vIdxX = round(linspace(1,vSz(2),iNumBlkX));
vIdxY = round(linspace(1,vSz(1),iNumBlkY));

% Fix vectors if only single block
if length(vIdxX) == 1
    vIdxX = [1 vIdxX];
end
if length(vIdxY) == 1
    vIdxY = [1 vIdxY];
end

% Initialize
if lMat
    mGrid = [];
else
    switch strClass
        case 'logical'
            mGrid = false(vSz);
        case 'uint8'
            mGrid = zeros(vSz);
        case 'double'
            mGrid = NaN(vSz);
        otherwise
            error('Invalid class type.');
    end
end

% Check whether user-specified class matches the matfile variable class
if lMat
    try
        strClassF = class(objM.(strField));
    catch
        strClassF = strClass;
    end
    if ~strcmp(strClass,strClassF)
        strClass = strClassF;
        warning('block_process:class_match',  ...
            'Specified class and matfile variable class do not match.')
    end
end

% Loop through each block
for iY = 1:length(vIdxY)-1
    for iX = 1:length(vIdxX)-1
        
        % Input for user-defined function
        sInput.indexX = vIdxX;
        sInput.indexY = vIdxY;
        sInput.countX = iX;
        sInput.countY = iY;
        
        % Call user-defined function for current block
        mBlock = hFun(sInput);
        
        % Spatial referencing for current block
        vIdxXb = vIdxX(iX):vIdxX(iX+1);
        vIdxYb = vIdxY(iY):vIdxY(iY+1);

        % Skip if user-defined function returned the wrong size
        if any(size(mBlock) ~= [length(vIdxYb) length(vIdxXb)])
            warning('MATLAB:wrongBlockSize', ...
                'Function returned the wrong block size.')
            continue
        end
        
        % Save block       
        if lMat           
            eval(['objM.' strField ...
                '(vStartIdx(1)+vIdxYb-1,vStartIdx(2)+vIdxXb-1)' ...
                ' = cast(mBlock,strClass);'])
        else
            mGrid(vIdxYb,vIdxXb) = cast(mBlock,strClass);
        end
        
    end
end