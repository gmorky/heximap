% Main script for Hexagon image stitching

% Waitbar handle
hW = waitbar(0/6,'please wait...');

% Make sure both halves are selected for all images
if ~iscell(cFile) || mod(numel(cFile),2) ~= 0
    error('Must select both "a" and "b" halves for each Hexagon image.')
end

% Detect image corners
cFileA = cFile(1:2:end); cFileB = cFile(2:2:end);
[cCornersA,cCornerFigsA] = cellfun(@(x) stiGetCorners(strPath,x,'a',hW), ...
    cFileA,'Uni',0);
[cCornersB,cCornerFigsB] = cellfun(@(x) stiGetCorners(strPath,x,'b',hW), ...
    cFileB,'Uni',0);

% Stitch together the image halves
cM = cellfun(@(a,b,c,d) stiStitch(strPath,a,b,c,d,hW), ...
    cFileA,cFileB,cCornersA,cCornersB,'Uni',0);

% Save lower resolution copies of the images
cellfun(@(x) stiResize(x,hW),cM,'Uni',0);

% Save additional image info
cellfun(@(a,b,c) stiSaveInfo(a,b,c),cM,cCornerFigsA,cCornerFigsB,'Uni',0);

% Update waitbar
try
waitbar(6/6,hW,'finished.')
catch
end
