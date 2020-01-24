function objM = stiStitch(strPath,strFileL,strFileR,mCornL,mCornR,hW)
% Detect and match feature points, stitch scanned Hexagon image halves

% Get image info
sInfoL = imfinfo(strcat([strPath strFileL]));
sInfoR = imfinfo(strcat([strPath strFileR]));

% Make sure correct structures are being accessed
if numel(sInfoL) > 1
    sInfoL = sInfoL(1);
end
if numel(sInfoR) > 1
    sInfoR = sInfoR(1);
end

% Create mat file to store images
[~,strN,~] = fileparts(sInfoL.Filename);
strFile = strcat([strPath strN(1:end-2)]);
if exist([strFile '.mat'],'file')
    delete([strFile '.mat'])
end
objM = matfile(strFile,'Writable',true);

% Update waitbar
try
waitbar(2/6,hW,'saving left half of image in mat file...')
catch
end

% Save left image half in the mat file
saveLeftImageHalf(sInfoL,mCornL,objM)

% Update waitbar
try
waitbar(3/6,hW,'detecting and matching ORB features...')
catch
end

% Estimate transformation for right image half to match left image half
mT = estimateTransform(sInfoR,mCornR,objM);

% Update waitbar
try
waitbar(4/6,hW,'saving right half of image in mat file...')
catch
end

% Save right image half in the mat file
saveRightImageHalf(sInfoR,mCornR,objM,mT)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = saveLeftImageHalf(sInfoL,mCornL,objM)
% Save left image half in the mat file

% Make index vectors
dNumSections = 10;
vR = mCornL(:,2)';
vRi = vR-mCornL(3)+1;
vC = round(linspace(mCornL(1),mCornL(2),dNumSections));

% Loop for each image section
for i = 1:length(vC)-1

    % Save image section in mat file
    vCi = [vC(i) vC(i+1)]-mCornL(1)+1;
    objM.Image(vRi(1):vRi(2),vCi(1):vCi(2)) = imread(sInfoL.Filename, ...
        'PixelRegion',{vR,[vC(i) vC(i+1)]});
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mT = estimateTransform(sInfoR,mCornR,objM)
    
% Get size of left image half
[iH,iW] = size(objM,'Image');

% Define pixel region boundaries
iX = round(iW/8);
vY = round(linspace(1,iH,8));
iNumWin = length(vY)-1;

% Initialize
mL = []; mR = [];

% Loop through each pixel region
for i = 1:iNumWin
    
    try
    
    % Construct classes
    clDetector = cv.FeatureDetector('ORB','MaxFeatures',10000);
    clExtractor = cv.DescriptorExtractor('ORB');
    clMatcher = cv.DescriptorMatcher('BruteForce-Hamming'); 
    
    % Read the left image, detect keypoints and descriptors, make
    % keypoint matrix
    mPixRegL = [iW-iX vY(i); iW vY(i+1)];
    mI = objM.Image(mPixRegL(3):mPixRegL(4),mPixRegL(1):mPixRegL(2));
    sKeyPoints = clDetector.detect(mI);
    sDescriptorsL = clExtractor.compute(mI,sKeyPoints);
    mKeyPoints = [sKeyPoints.pt];
    mKeyPointsL = [mKeyPoints(1:2:end); mKeyPoints(2:2:end)];
    
    % Read the right image, detect keypoints and descriptors, make
    % keypoint matrix
    cPixRegR = {[vY(i) vY(i+1)]+mCornR(3) [1 iX]+mCornR(1)};
    mI = imread(sInfoR.Filename,'PixelRegion',cPixRegR);
    sKeyPoints = clDetector.detect(mI);
    sDescriptorsR = clExtractor.compute(mI,sKeyPoints);
    mKeyPoints = [sKeyPoints.pt];
    mKeyPointsR = [mKeyPoints(1:2:end); mKeyPoints(2:2:end)];
    clear mI sKeyPoints mKeyPoints

    % Match the keypoints
    sMatches = clMatcher.match(sDescriptorsL,sDescriptorsR);

    % Make vectors of matched points
    vMatchIdxL = [sMatches.queryIdx]+1;
    vMatchIdxR = [sMatches.trainIdx]+1;
    mPtsL = mKeyPointsL(:,vMatchIdxL);
    mPtsR = mKeyPointsR(:,vMatchIdxR);
    
    % Use ransac to estimate inliers
    sRansac.numIter = 2000;
    sRansac.inlierDist = 1;
    [~,lIn] = estimateTransformRansac([mPtsL;ones(1,size(mPtsL,2))], ...
                                      [mPtsR;ones(1,size(mPtsR,2))], ...
                                      sRansac);
    mPtsL = mPtsL(:,lIn);
    mPtsR = mPtsR(:,lIn);
    
    % Convert points to full image coordinates
    mPtsL(1,:) = mPtsL(1,:) + mPixRegL(1) - 1;
    mPtsL(2,:) = mPtsL(2,:) + mPixRegL(3) - 1;
    mPtsR(1,:) = mPtsR(1,:) + cPixRegR{2}(1) - 1;
    mPtsR(2,:) = mPtsR(2,:) + cPixRegR{1}(1) - 1;
    
    % Make points relative to axis with origin at lower left corner
    mPtsL(2,:) = iH - mPtsL(2,:) + 1;
    mPtsR(2,:) = sInfoR.Height - mPtsR(2,:) + 1;
    
    % Save points
    mL = [mL mPtsL];
    mR = [mR mPtsR]; 
    
    catch
    end

end

% Estimate geometric transform
sRansac.NumIter = 2000;
sRansac.InlierDist = 1;
[mT,~] = estimateTransformRansac([mL;ones(1,size(mL,2))], ...
                                 [mR;ones(1,size(mR,2))], ...
                                 sRansac);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = saveRightImageHalf(sInfoR,mCornR,objM,mT)
    
% Get image size
[iH,iW] = size(objM,'Image');

% Spatial referencing info for image to be transformed (right half)
sR.DeltaLon = 1;
sR.DeltaLat = -1;
sR.Lonlim = [1 sInfoR.Width] + [-0.5 0.5];
sR.Latlim = [1 sInfoR.Height] + [-0.5 0.5];
sInfoR.SpatialRef = sR;

% Spatial referencing vectors for saving in mat file
mWin = round(mT * [mCornR(2,:) 1 1; mCornR(2) 1 1 1]');
vX = iW:max(mWin(1,:));
vY = fliplr(1:iH);

% Transform the right image half and save it in the mat file
sParams.blockSize = 1000;  
sParams.matSource = objM.Properties.Source;
sParams.matField = 'Image';  
sParams.matStartIndex = [min(vY) min(vX)];
sParams.transform = mT;
grid2grid(sInfoR,vX,vY,sParams);

end
end
