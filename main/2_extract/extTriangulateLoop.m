function [] = extTriangulateLoop(strHexPath,hW)
% Triangulate points for each region

% Update waitbar
try
waitbar(7/8,hW,'triangulating points...')
catch
end

% Initialize
cFileL = cellfun(@(x) matfile(x,'Writable',true), ...
    getFiles(strHexPath,'Left.mat'),'Uni',0);
cFileR = cellfun(@(x) matfile(x,'Writable',true), ...
    getFiles(strHexPath,'Right.mat'),'Uni',0);
vROI = cell2mat(cellfun(@(x) x.RegionID,cFileL,'Uni',0));
vUniqueROI = unique(vROI);
iNumROI = length(vUniqueROI);

% Loop through each ROI
for iR = 1:iNumROI
    
    try
    
    % Initialize
    cReg = {num2str(iR) num2str(iNumROI)};
    
    % Get windows belonging to current ROI    
    [cL,cR] = extGetROI(cFileL,cFileR,vROI,vUniqueROI(iR));

    % Loop through windows belonging to current ROI
    for iW = 1:numel(cL)
        
        % Update command window
        disp(['triangulating points for region ' ...
            num2str(vUniqueROI(iR)) ' window ' num2str(iW) '...'])
        
        % Triangulate the points
        extTriangulate(cL{iW},cR{iW},hW,cReg);
        
    end
    
    catch objExc

    warning(objExc.message)
    warning(['An error occurred during triangulation. ' ...
        'Skipping region ' num2str(vUniqueROI(iR)) '...'])
    for iW = 1:numel(cL)
        cL{iW}.Error = objExc.message;
        cR{iW}.Error = objExc.message;
    end
    
    end
end
