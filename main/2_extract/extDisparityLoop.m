function [] = extDisparityLoop(cM,strHexPath,cWindow,strRes,iBlkSz,hW)
% Rectify stereo images and compute disparity maps for each window

% Update waitbar
try
waitbar(5/8,hW,'rectifying images and computing disparity maps...')
catch
end

% Warn user if files exist from a previous run
cFileL = getFiles(strHexPath,'Left.mat');
cFileR = getFiles(strHexPath,'Right.mat');
if ~isempty(cFileL) || ~isempty(cFileR)
    strQ = questdlg([ ...
        'Existing ''Left.mat'' and ''Right.mat'' files ' ...
        'in folder will cause errors. ' ...
        'Please move them, then click continue.'], ...
        '','Continue','Cancel','Cancel');                 
    if ~strcmpi(strQ,'Continue')
        error('Process cancelled by user.')
    end
end

% Initialize
iNumWin = numel(cWindow);

% Loop through each window
for iW = 1:iNumWin
        
    try
    
    % Update command window
    disp(['processing window ' num2str(iW) '...'])
    
    % Delete old directory and make new one for saving data
    strSavePath = [strHexPath num2str(iW) '\'];
    if exist(strSavePath,'dir')
        warning('off','MATLAB:RMDIR:RemovedFromPath');
        rmdir(strSavePath,'s');
        warning('on','MATLAB:RMDIR:RemovedFromPath');
    end
    mkdir(strHexPath,num2str(iW))
    
    % Get full Hexagon image mat files
    objML = cM{1};objMR = cM{2};
  
    % Make mat file for left Hexagon window
    strFile = strcat([strSavePath 'Left' '.mat']);
    objL = matfile(strFile,'Writable',true);
    
    % Make mat file for right Hexagon window
    strFile = strcat([strSavePath 'Right' '.mat']);
    objR = matfile(strFile,'Writable',true);
    
    % Save camera matrices in new mat files
    objL.IntrinsicMatrix = objML.IntrinsicMatrix;
    objR.IntrinsicMatrix = objML.IntrinsicMatrix;
    objL.PoseMatrix = objML.LeftPoseMatrix;
    objR.PoseMatrix = objML.RightPoseMatrix;
    
    % Save source image info in mat files
    objL.SourceImageInfo = struct('Path',objML.Properties.Source, ...
        'SpatialTrans',objML.SpatialTrans);
    objL.SourceImage = objML.Image10;
    objR.SourceImageInfo = struct('Path',objMR.Properties.Source, ...
        'SpatialTrans',objMR.SpatialTrans);
    objR.SourceImage = objMR.Image10;
    
    % Save window and ROI info in mat files
    objL.Window = cWindow{iW}.left;
    objR.Window = cWindow{iW}.right;
    objL.WindowID = iW;
    objR.WindowID = iW;
    objL.RegionID = cWindow{iW}.region;
    objR.RegionID = cWindow{iW}.region;
    
    % Initialize
    cWin = {num2str(iW) num2str(iNumWin)};
    objL.Accuracy = struct();
                  
    % Read Hexagon images    
    extReadImage(objML,objL,hW,cWin);
    extReadImage(objMR,objR,hW,cWin);

    % Rectify the stereo images
    extStereoRect(objL,objR,strSavePath,hW,cWin);
    
    % Compute disparity map 
    extDisparity(objL,objR,strRes,iBlkSz,hW,cWin);

    catch objExc

    warning(objExc.message)
    warning(['An error occurred during stereo rectification or ' ...
        'disparity map computation. Skipping window ' num2str(iW) '...'])
    objL.Error = objExc.message;
    objR.Error = objExc.message;
    
    end
end
