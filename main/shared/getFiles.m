function cFiles = getFiles(strDir,strExt)

% Get data for the current directory
sDirData = dir(strDir);

% Get index for directories
lDirIdx = [sDirData.isdir];

% Get index for files with correct extension
cDir = {sDirData.name}';
lFileIdx = false(size(cDir));
for i = 1:length(cDir)
    if strfind(cDir{i},strExt) == length(cDir{i})-length(strExt)+1
        lFileIdx(i) = true;
    end
end

% Get a list of the files
cFiles = {sDirData(lFileIdx).name}';

% Prepend path to files
if ~isempty(cFiles)
    cFiles = cellfun(@(x) fullfile(strDir,x),cFiles,'Uni',0);
end

% Get a list of the subdirectories
cSubDirs = {sDirData(lDirIdx).name};

% Find index of subdirectories that are not '.' or '..'
lValidIdx = ~ismember(cSubDirs,{'.','..'});

% Loop over valid subdirectories
for iDir = find(lValidIdx)
    
    % Get the subdirectory path
    strNextDir = fullfile(strDir,cSubDirs{iDir});
    
    % Recursively call getFiles
    cFiles = [cFiles; getFiles(strNextDir,strExt)];
    
end
