function lGrid = polygons2grid(polygons,vX,vY,varargin)

% Size of output grid
vSz = [length(vY) length(vX)];

% Default parameters
if ~isempty(varargin)
    sParams = varargin{1};
else
    sParams.blockSize = max(vSz)+1;
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
    sParams.blockSize = max(vSz)+1;    
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

% Initialize logical array in mat file
if lMat
    objM = matfile(sMatSaveInfo.file,'Writable',true);
    objM.(sMatSaveInfo.field) = false(2,2);
    clear objM
end

if ischar(polygons)

    % Get paths to shapefiles
    cFiles = getFiles(polygons,'.shp');
    cFiles = cFiles(cellfun(@isempty,strfind(cFiles,'.xml')));

else
    
    % Assign shapefile
    cFiles{1} = polygons;
    
end

% Function handle for processing blocks
hFun = @(x) makeBlock(x,vX,vY,cFiles);

% Call block processing function
warning('off','MATLAB:wrongBlockSize');
if lMat
    lGrid = blockProcess(vSz,sParams.blockSize,'logical',hFun,sMatSaveInfo);
else
    lGrid = blockProcess(vSz,sParams.blockSize,'logical',hFun);
end
warning('on','MATLAB:wrongBlockSize');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lBlock = makeBlock(sInput,vX,vY,cFiles)
    
% Get input from block processing function
vIdxX = sInput.indexX;
vIdxY = sInput.indexY;
iX = sInput.countX;
iY = sInput.countY;

% Spatial referencing grids for new image block
vBlkX = vX(vIdxX(iX):vIdxX(iX+1));
vBlkY = vY(vIdxY(iY):vIdxY(iY+1));

% Initialize
lBlock = false(length(vBlkY),length(vBlkX));

% Loop for each file
for i = 1:numel(cFiles)
    
    % Read shapefile if user supplied the path
    if ischar(cFiles{i}) 
        sP = shaperead(cFiles{i},'BoundingBox', ...
            [vBlkX(1) vBlkX(end); vBlkY(end) vBlkY(1)]','Attributes',{}); 
    else
        sP = cFiles{i};
    end
    
    % Make block for this shapefile
    lP = makeGrid(sP,vBlkX,vBlkY);

    % Update block
    lBlock(lP) = true;
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lP = makeGrid(sM,vX,vY)

% Initialize mask
lP = false([length(vY) length(vX)]);

% Loop for each shapefile polygon
for iP = 1:length(sM)  
    
    % Make Nx2 matrix of outer polygon vertices
    vIdx = find(isnan(sM(iP).X)|isnan(sM(iP).Y));
    vXo = sM(iP).X(1:vIdx(1)-1);
    vYo = sM(iP).Y(1:vIdx(1)-1);
    mOuter = [vXo(:) vYo(:)];
    
    % Make outer polygon mask
    lIn = roipoly(vX,vY,lP,mOuter(:,1),mOuter(:,2));
    
    % Loop for each hole in polygon
    for iH = 1:length(vIdx)-1
        
        % Make Nx2 matrix of inner polygon (hole) vertices
        mInner = [sM(iP).X(vIdx(iH)+1:vIdx(iH+1)-1)' ...
                  sM(iP).Y(vIdx(iH)+1:vIdx(iH+1)-1)'];
        
        % Unmask pixels within polygon holes
        lHoles = roipoly(vX,vY,lP,mInner(:,1),mInner(:,2));     
        lIn(lHoles & lIn) = false;
        
    end
    
    % Save in polygon matrix
    lP(lIn) = true;

end
end
end