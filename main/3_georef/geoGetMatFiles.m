function [cL,cR] = geoGetMatFiles(strWinPath,hW)
% Get mat files containing image windows

% Update waitbar
try
waitbar(2/8,hW,'checking .mat files...')
catch
end

% Get paths to mat files
cL = cellfun(@(x) matfile(x,'Writable',true),getFiles(strWinPath, ...
    'Left.mat'),'Uni',0);
cR = cellfun(@(x) matfile(x,'Writable',true),getFiles(strWinPath, ...
    'Right.mat'),'Uni',0);

% Make sure corresponding left and right images exist for each window
if numel(cL) ~= numel(cR)
    error('Number of left and right .mat files is not equal.')
end

% Remove any windows where extraction was not completed
lIn = true(size(cL));
for i = 1:numel(cL)
    if isempty(whos(cL{i},'TriangulatedPoints'))
        lIn(i) = false;
        warning(['Extraction was not completed for window ' ...
            num2str(cL{i}.WindowID) '. Skipping.'])
    end
end
cL = cL(lIn);
cR = cR(lIn);

% Make sure at least one valid window exists
if isempty(cL)
    error('No valid windows found in the specified folder.')
end

% Sort the windows
[~,vIdx] = sort(cell2mat(cellfun(@(x) x.WindowID,cL,'Uni',0)));
cL = cL(vIdx); cR = cR(vIdx);
