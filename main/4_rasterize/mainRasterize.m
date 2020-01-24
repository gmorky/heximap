% Main script for exporting raster geotiff files

% Waitbar handle
hW = waitbar(0/1,'please wait...');

% Get mat files containing data for image windows
cL = rasGetMatFiles(strWinPath,hW);

% Make directories for saving rasters
mkdir(strWinPath,'dems')
mkdir(strWinPath,'images')

% Loop through each window
for iW = 1:numel(cL)
    
    try
    
    % Initialize
    cWin = {num2str(iW) num2str(numel(cL))};
    
    % Rasterize the DEM
    rasDem(cL{iW},strWinPath,hW,cWin);
    
    % Rasterize the orthoimage
    rasOrtho(cL{iW},strWinPath,hW,cWin);
    
    catch objExc

    warning(objExc.message)
    warning(['An error occurred. Skipping window ' num2str(iW) '...'])
    cL{iW}.Error = objExc.message;
    
    end    
end

% Update waitbar
try
waitbar(1/1,hW,'finished.')
catch
end
