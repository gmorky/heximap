function [] = extDisparity(objL,objR,strRes,hW,cWin)
% Compute disparity maps

% Update waitbar
try
set(get(findobj(hW,'type','axes'),'title'), 'string', ...
    ['window ' cWin{1} ' of ' cWin{2} ': computing disparity map...'])
pause(0.1)
catch
end

% Filter images
sOpt.histmatch = true;
sOpt.adapthisteq = true;
sOpt.wiener2 = false;
[objL.RectImageFilt,objR.RectImageFilt] = extFilterImages( ...
    objL.RectImage,objR.RectImage,0,sOpt);

% Spatial resolution of output (disparity map resampling factors)
switch strRes
    case 'Very High'
        vScale = [16 8 4 2 1];
    case 'High'
        vScale = [16 8 4 2];
    case 'Medium'
        vScale = [16 8 4];
    case 'Low'
        vScale = [16 8];
    otherwise
        error('Unknown scale tag.')
end

% Define paramaters for disparity algorithm
vD = [-8 8] * 30;
iBlkSz = 15;
sStereo.MinDisparity = vD(1);
sStereo.NumDisparities = vD(2)-vD(1);
sStereo.BlockSize = iBlkSz;
sStereo.P1 = 0;
sStereo.P2 = 32 * iBlkSz^2;
sStereo.Disp12MaxDiff = 1;
sStereo.PreFilterCap = 0;
sStereo.UniquenessRatio = 0;
sStereo.SpeckleWindowSize = 300;
sStereo.SpeckleRange = 1;
sStereo.Mode = 'SGBM';

% Outlier threshold
dThresh = 50;

% Loop through each scale
for i = 1:length(vScale)

    % Compute disparity map
    [mD,vScaleS] = sgbm(objL,objR,sStereo,vScale(i));
    
    % Multiply disparity map by scale factor
    mD = mD * vScaleS(2);    
    
    if i > 1
        
        % Remove outliers
        mLoRes = imresize(mLoRes,size(mD));
        mD(abs(mD-mLoRes)>dThresh) = NaN;
        
        if i < length(vScale)
            
            % Fill holes with lower resolution disparity map
            lFill = ~isnan(mLoRes) & isnan(mD);
            mD(lFill) = mLoRes(lFill);
            
        end
    end
    
    if i < length(vScale)
        
        % Save lower resolution disparity map for next iteration
        mLoRes = mD;   
        
    else
        
        % Find non-empty disparity pixels, then save disparity map
        lIn = ~isnan(mD(:));
        objL.Disparity = mD;
        objL.DisparityScale = vScaleS;
        
    end
    
    % Clear to conserve memory
    clear mD

end

% Get data
mHexWinL = objL.Window;
mHexWinR = objR.Window;
mProjWin = objL.RectWindow;
mH1 = objL.Homography;
mH2 = objR.Homography;

% Make Nx3 matrix of pre-rectified left image points
[mX,mY] = meshgrid(mProjWin(1):mProjWin(2),mProjWin(3):mProjWin(4));
mX = imresize(mX,size(objL,'RectImageFilt')./vScaleS);
mY = imresize(mY,size(objL,'RectImageFilt')./vScaleS);
mPts = mH1 \ [mX(lIn) mY(lIn) ones(sum(lIn),1)]';
mPts(1,:) = mPts(1,:) ./ mPts(3,:) + mHexWinL(1) - 1;
mPts(2,:) = mPts(2,:) ./ mPts(3,:) + mHexWinL(3) - 1;
mPts(3,:) = mPts(3,:) ./ mPts(3,:);
objL.ImagePoints = mPts;
clear mPts

% Make Nx3 matrix of pre-rectified right image points
[mX,mY] = meshgrid(mProjWin(1):mProjWin(2),mProjWin(3):mProjWin(4));
mX = imresize(mX,size(objR,'RectImageFilt')./vScaleS);
mY = imresize(mY,size(objR,'RectImageFilt')./vScaleS);
mD = objL.Disparity;
mPts = mH2 \ [mX(lIn)-mD(lIn) mY(lIn) ones(sum(lIn),1)]';
clear mD
mPts(1,:) = mPts(1,:) ./ mPts(3,:) + mHexWinR(1) - 1;
mPts(2,:) = mPts(2,:) ./ mPts(3,:) + mHexWinR(3) - 1;
mPts(3,:) = mPts(3,:) ./ mPts(3,:);
objR.ImagePoints = mPts;
clear mPts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mD,vScaleS] = sgbm(objL,objR,sStereo,dScale)
% Compute disparity map using semi-global block matching

% Determine closest integer size of image at given scale
vSz = size(objL,'RectImageFilt');
vNewSz = round(vSz/dScale);

% Determine vertical and horizontal scale factors
vScaleS = vSz./vNewSz;

% Set disparity search range for current scale
vDs(1) = sStereo.MinDisparity;
vDs(2) = vDs(1) + sStereo.NumDisparities;
vDs = round(vDs/dScale/8)*8;
if any(abs(vDs) < 8)
    vDs = [-8 8];
end
    
% Set parameters
sStereo = cv.StereoSGBM('MinDisparity', vDs(1), ...
                        'NumDisparities',vDs(2)-vDs(1), ...
                        'BlockSize',sStereo.BlockSize, ...
                        'P1',sStereo.P2, ...
                        'P2',sStereo.P1, ...
                        'Disp12MaxDiff',sStereo.Disp12MaxDiff, ...
                        'PreFilterCap',sStereo.PreFilterCap, ...
                        'UniquenessRatio',sStereo.UniquenessRatio, ...
                        'SpeckleWindowSize',sStereo.SpeckleWindowSize, ...
                        'SpeckleRange',sStereo.SpeckleRange, ...
                        'Mode',sStereo.Mode); 

% Compute disparity map
mD = double(sStereo.compute(imresize(objL.RectImageFilt,vNewSz), ...
    imresize(objR.RectImageFilt,vNewSz))) / 16;

% Remove invalid disparities
mD(mD < vDs(1) | mD > vDs(2)) = NaN;

% Remove erroneous regions around edges of disparity image
lMask = imresize(objL.RectImage,vNewSz) == 0 | ...
        imresize(objR.RectImage,vNewSz) == 0;
lMask(imclearborder(lMask)) = false;
lMask = imdilate(lMask,ones(ceil(100/dScale)));
mD(lMask) = NaN;

end
end