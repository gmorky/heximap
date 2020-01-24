function [] = extBundleAdjustLoop(strHexPath,hW)
% Perform bundle adjustment for region

% Update waitbar
try
waitbar(6/8,hW,'performing bundle adjustment...')
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
    
    % Update command window
    disp(['performing bundle adjustment for region ' ... 
        num2str(vUniqueROI(iR)) '...'])
        
    % Initialize
    cReg = {num2str(iR) num2str(iNumROI)};
    
    % Get windows belonging to current ROI    
    [cL,cR] = extGetROI(cFileL,cFileR,vROI,vUniqueROI(iR));
    
    % Perform bundle adjustment by minimizing reprojection error
    extBundleAdjust(strHexPath,cL,cR,hW,cReg);
    
    catch objExc

    warning(objExc.message)
    warning(['An error occurred during bundle adjustment. ' ...
        'Skipping region ' num2str(vUniqueROI(iR)) '...'])
    for iW = 1:numel(cL)
        cL{iW}.Error = objExc.message;
        cR{iW}.Error = objExc.message;
    end
    
    end   
end

% Loop through each ROI
for iR = 1:iNumROI
    
    try
    
    % Initialize
    cReg = {num2str(iR) num2str(iNumROI)};
    
    % Get windows belonging to current ROI    
    [cL,cR] = extGetROI(cFileL,cFileR,vROI,vUniqueROI(iR));
    
    % Try bundle adjustment again for any previous inaccurate windows. This
    % time, can use better starting guess for solver (using solutions from
    % other more accurate windows)
    sA = cL{1}.Accuracy;
    if sA.BundleAdjust > 10000
        
        % Update command window
        disp(['redoing bundle adjustment for region ' ...
            num2str(vUniqueROI(iR)) '...'])
        
        % Redo bundle adjustment
        extBundleAdjust(strHexPath,cL,cR,hW,cReg);
        
    end
    
    catch objExc

    warning(objExc.message)
    warning(['An error occurred during bundle adjustment redo. ' ...
        'Skipping region ' num2str(vUniqueROI(iR)) '...'])
    for iW = 1:numel(cL)
        cL{iW}.Error = objExc.message;
        cR{iW}.Error = objExc.message;
    end
    
    end  
    
end
