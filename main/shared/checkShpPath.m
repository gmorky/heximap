function [] = checkShpPath(strPath)
% Check given folder containing shapefiles for errors

% Get paths to shapefiles
cFiles = getFiles(strPath,'.shp');

% Make sure xml files are not included
try
    cFiles = cFiles(cellfun(@isempty,strfind(cFiles,'.xml')));
catch
end

% Throw error is no valid shapefiles are found in the folder
if isempty(cFiles)
    error('No valid shapefiles were found in the specified folder.')
end

% Loop through and check each shapefile
for i = 1:numel(cFiles)
    checkInput(cFiles{i});
end
