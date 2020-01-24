function [] = extStereoRect(objL,objR,strSavePath,hW,cWin)
% Rectify stereo images so corresponding features line up along rows. This
% reduces the disparity search to one dimension. To get a better
% distribution of point matches across the images, each image is divided
% into blocks and ORB points are detected seperately for each block.

% Update waitbar
try
set(get(findobj(hW,'type','axes'),'title'), 'string', ...
    ['window ' cWin{1} ' of ' cWin{2} ': computing epipolar images...'])
pause(0.1)
catch
end

% Filter the images
sOpt.histmatch = true;
sOpt.adapthisteq = true;
sOpt.wiener2 = true;
[objL.ImageFilt,objR.ImageFilt] = extFilterImages( ...
    objL.Image,objR.Image,0,sOpt);

% Parameters
mPtsL = []; mPtsR = [];
sParams.maxFeatures = 10000;
sParams.numKeep = 100;
sParams.maxError = 6;
dBlkSz = round(mean(size(objL,'Image')/4));

% Process image blocks to get initial matches
warning('off','MATLAB:wrongBlockSize');
blockProcess(size(objL,'Image'),dBlkSz,'double', ...
    @(x) block1(x,objL,objR,sParams));

% Compute initial homographies
[mH1,mH2] = homographies(size(objL,'Image'),sParams.maxError);

% Parameters
mPtsL = []; mPtsR = [];
sParams.maxError = 3;

% Process image blocks get improved distribution of matches
blockProcess(size(objL,'Image'),dBlkSz,'double', ...
    @(x) block2(x,objL,objR,mH1,mH2));
warning('on','MATLAB:wrongBlockSize');

% Compute new homographies
[mH1,mH2] = homographies(size(objL,'Image'),sParams.maxError);

% Transform images
transform(objL,objR,strSavePath,mH1,mH2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mBlock = block1(sInput,objL,objR,sParams)

% Null output
mBlock = [];
    
% Get input from block processing function
vIdxX = sInput.indexX;
vIdxY = sInput.indexY;
iX = sInput.countX;
iY = sInput.countY;

% Image size
vSz = size(objL,'Image');

% Read left image block
mBlkL = [vIdxX(iX) vIdxY(iY); vIdxX(iX+1) vIdxY(iY+1)];
mI1 = objL.ImageFilt(mBlkL(3):mBlkL(4),mBlkL(1):mBlkL(2));

% Read right image block
iB = 200;
mBlkR = mBlkL + [-iB -iB; iB iB];
mBlkR(mBlkR < 1) = 1;
mBlkR(mBlkR(:,1) > vSz(2),1) = vSz(2);
mBlkR(mBlkR(:,2) > vSz(1),2) = vSz(1);
mI2 = objR.ImageFilt(mBlkR(3):mBlkR(4),mBlkR(1):mBlkR(2));

% Detect keypoints and compute descriptors for left image block
sKeypointsL = keypoints(mI1,sParams.maxFeatures);
mDescriptorsL = descriptors(mI1,sKeypointsL);

% Detect keypoints and compute descriptors for right image block
sKeypointsR = keypoints(mI2,sParams.maxFeatures);
mDescriptorsR = descriptors(mI2,sKeypointsR);

% Find matching keypoints
[mPtsLb,mPtsRb] = matches(sKeypointsL,sKeypointsR, ...
    mDescriptorsL,mDescriptorsR,sParams.numKeep);

% Convert left image block keypoints to full window coordinates
mPtsLb(1,:) = mPtsLb(1,:) + mBlkL(1) - 1;
mPtsLb(2,:) = mPtsLb(2,:) + mBlkL(3) - 1;

% Convert right image block keypoints to full window coordinates
mPtsRb(1,:) = mPtsRb(1,:) + mBlkR(1) - 1;
mPtsRb(2,:) = mPtsRb(2,:) + mBlkR(3) - 1;

% Save output
mPtsL = [mPtsL mPtsLb];
mPtsR = [mPtsR mPtsRb];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sKeypoints = keypoints(mI,iMaxFeatures)
% Detect keypoints
    
% Construct classes
clDetector = cv.FeatureDetector('ORB','MaxFeatures',iMaxFeatures);

% Detect keypoints
sKeypoints = clDetector.detect(mI);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mDescriptors = descriptors(mI,sKeypoints)
% Extract descriptors
        
% Construct classes
clExtractor = cv.DescriptorExtractor('ORB');

% Extract descriptors
mDescriptors = clExtractor.compute(mI,sKeypoints);
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sKeypointsL,sKeypointsR] = matches(sKeypointsL,sKeypointsR, ...
        mDescriptorsL,mDescriptorsR,iNumKeep)
% Find matching keypoints

% Construct class
clMatcher = cv.DescriptorMatcher('BruteForce-Hamming');

% Match the keypoints
sMatches = clMatcher.match(mDescriptorsL,mDescriptorsR);    
[~,vIdx] = sort([sMatches.distance],'ascend');

try

% Only keep the strongest matches
sMatches = sMatches(vIdx(1:iNumKeep));

catch
end

try

% Make vectors of matched points
vMatchIdxL = [sMatches.queryIdx]+1;
vMatchIdxR = [sMatches.trainIdx]+1;
sKeypointsL = sKeypointsL(:,vMatchIdxL)';
sKeypointsR = sKeypointsR(:,vMatchIdxR)';

% Output
sKeypointsL = [sKeypointsL.pt];
sKeypointsL = [sKeypointsL(1:2:end); sKeypointsL(2:2:end)];
sKeypointsR = [sKeypointsR.pt];
sKeypointsR = [sKeypointsR(1:2:end); sKeypointsR(2:2:end)];

catch
    
% Empty output
sKeypointsL = [];
sKeypointsR = [];
   
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mH1,mH2] = homographies(vSz,iMaxError)

% Initialize
iNumAttempts = 100;
iNumPts = size(mPtsL,2);
vError = zeros(1,iNumAttempts);
cH1 = cell(1,iNumAttempts);
cH2 = cell(1,iNumAttempts);
cPtsL = cell(1,iNumAttempts);
cPtsR = cell(1,iNumAttempts);

% Attempts
for i = 1:iNumAttempts
    
    % Choose random subset of points
    vIdx = randperm(iNumPts,round(iNumPts*0.95));
    mPtsLp = mPtsL(:,vIdx);
    mPtsRp = mPtsR(:,vIdx);
        
    % Remove outliers using ransac
    [~,lIn] = cv.findFundamentalMat( ...
        num2cell(mPtsLp,1), ...
        num2cell(mPtsRp,1), ...
        'Method','Ransac','RansacReprojThreshold',iMaxError,'Confidence',0.99);
    mPtsLp = mPtsLp(:,logical(lIn));
    mPtsRp = mPtsRp(:,logical(lIn));
    
    % Compute fudamental matrix using 8-point algorithm
    [mF,~] = cv.findFundamentalMat( ...
        num2cell(mPtsLp,1), ...
        num2cell(mPtsRp,1), ...
        'Method','8Point');

    % Compute homographies
    [mH1,mH2] = cv.stereoRectifyUncalibrated( ...
        num2cell(mPtsLp,1), ...
        num2cell(mPtsRp,1), ...
        mF,vSz);
    
    % Transform points and compute number of points with vertical error
    % less than the specified max
    vError(i) = sum(error(mH1,mH2) < iMaxError);
    
    % Store data in cell arrays
    cH1{i} = mH1;
    cH2{i} = mH2;
    cPtsL{i} = mPtsLp;
    cPtsR{i} = mPtsRp;

end

% Choose best solution
[~,iIdx] = max(vError);
mH1 = cH1{iIdx};
mH2 = cH2{iIdx};
mPtsL = cPtsL{iIdx};
mPtsR = cPtsR{iIdx};

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vError = error(mH1,mH2)

% Transform point matches using homographies
mPtsLp = mH1 * [mPtsL; ones(1,size(mPtsL,2))];
mPtsRp = mH2 * [mPtsR; ones(1,size(mPtsR,2))];
mPtsLp = mPtsLp ./ repmat(mPtsLp(3,:),[3 1]);
mPtsRp = mPtsRp ./ repmat(mPtsRp(3,:),[3 1]);

% Vertical error
vError = mPtsLp(2,:) - mPtsRp(2,:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = transform(objL,objR,strSavePath,mH1,mH2)
    
% Image size
vSz = size(objL,'Image');
        
% Project corners for left image window
mProjWin1 = mH1 * [1 1 1; vSz(2) 1 1; vSz(2) vSz(1) 1; 1 vSz(1) 1]';
mProjWin1 = (mProjWin1(1:2,:) ./ repmat(mProjWin1(3,:),[2 1]))';

% Project corners for right image window
mProjWin2 = mH2 * [1 1 1; vSz(2) 1 1; vSz(2) vSz(1) 1; 1 vSz(1) 1]';
mProjWin2 = (mProjWin2(1:2,:) ./ repmat(mProjWin2(3,:),[2 1]))';

% Make final (combined) window
mProjWin = round([min([mProjWin1;mProjWin2]);max([mProjWin1;mProjWin2])]);
vSzP = diff(mProjWin)+1;

% Make sure transformed image will be a reasonable size
if any(vSzP > vSz*1.5)
    error(['Rectified image is too large, (' num2str(vSzP(1)) ' by ' ...
        num2str(vSzP(2))  ' pixels), ' ...
        'indicating poor stereo rectification results.'])
end

% Make grids and transform left image
[mPts,mY] = meshgrid(mProjWin(1):mProjWin(2),mProjWin(3):mProjWin(4));
mPts = mH1 \ [mPts(:) mY(:) ones(prod(vSzP),1)]'; clear mY
mPts(1,:) = mPts(1,:) ./ mPts(3,:);
mPts(2,:) = mPts(2,:) ./ mPts(3,:);
mPts = mPts(1:2,:);
objL.RectImage = cv.remap(objL.Image, ...
    reshape(mPts(1,:),fliplr(vSzP)),reshape(mPts(2,:),fliplr(vSzP)));

% Make grids and transform right image
[mPts,mY] = meshgrid(mProjWin(1):mProjWin(2),mProjWin(3):mProjWin(4));
mPts = mH2 \ [mPts(:) mY(:) ones(prod(vSzP),1)]'; clear mY
mPts(1,:) = mPts(1,:) ./ mPts(3,:);
mPts(2,:) = mPts(2,:) ./ mPts(3,:);
mPts = mPts(1:2,:);
objR.RectImage = cv.remap(objR.Image, ...
    reshape(mPts(1,:),fliplr(vSzP)),reshape(mPts(2,:),fliplr(vSzP)));
clear mPts

% Transform points and compute vertical error
vError = error(mH1,mH2);

% Save point matches in mat files
objL.PointMatches = mPtsL;
objR.PointMatches = mPtsR;

% Save transformations info
objL.Homography = mH1;
objR.Homography = mH2;
objL.RectWindow = mProjWin;
objR.RectWindow = mProjWin;

% Save point match accuracies (vertical error) in mat file
sAccuracy = objL.Accuracy;
sAccuracy.PointMatches = vError;
objL.Accuracy = sAccuracy;

% Save left image ORB points
f = figure;
imagesc(objL.Image),colormap(bone(256)),hold on
scatter(mPtsL(1,:),mPtsL(2,:),abs(vError*100),[1 0.5 0])
title('Left Image ORB Matches')
axis equal, axis off
saveas(f,[strSavePath 'LeftMatches'],'fig')
close(f)   

% Save right image ORB points
f = figure;
imagesc(objR.Image),colormap(bone(256)),hold on
scatter(mPtsR(1,:),mPtsR(2,:),abs(vError*100),[1 0.5 0])
title('Right Image ORB Matches')
axis equal, axis off
saveas(f,[strSavePath 'RightMatches'],'fig')
close(f)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mBlock = block2(sInput,objL,objR,mH1,mH2)

% Null output
mBlock = [];
    
% Get input from block processing function
vIdxX = sInput.indexX;
vIdxY = sInput.indexY;
iX = sInput.countX;
iY = sInput.countY;

% Choose 25 random points from left image block
mPtsLb = [randi([vIdxX(iX) vIdxX(iX+1)],[1 25]); ...
          randi([vIdxY(iY) vIdxY(iY+1)],[1 25])];

% Transform keypoints to right image window coordinates
iNumPtsL = size(mPtsLb,2);
mPtsRb = mH2 \ mH1 * [mPtsLb;ones(1,iNumPtsL)];
mPtsRb = round(mPtsRb(1:2,:)./ repmat(mPtsRb(3,:),[2 1]));

% Read image windows
mI1 = objL.ImageFilt;
mI2 = objR.ImageFilt;

% Initialize
[iH,iW] = size(mI1);
vSearch = [200 20];

% Loop through each keypoint
for i = 1:iNumPtsL
    
    % Define sub-block indices for current keypoint
    mBlkL = round([mPtsLb(:,i)'-vSearch;mPtsLb(:,i)'+vSearch]);
    mBlkR = round([mPtsRb(:,i)'-vSearch;mPtsRb(:,i)'+vSearch]);
    
    % Skip if left sub-block is too close to edges
    if mBlkL(1) <  1 || mBlkL(3) <  1 || ...
       mBlkL(2) > iW || mBlkL(4) > iH
        continue
    end
    
    % Skip if right sub-block is too close to edges
    if mBlkR(1) <  1 || mBlkR(3) <  1 || ...
       mBlkR(2) > iW || mBlkR(4) > iH
        continue
    end

    % Get image sub-blocks
    mI1b = objL.ImageFilt(mBlkL(3):mBlkL(4),mBlkL(1):mBlkL(2));
    mI2b = objR.ImageFilt(mBlkR(3):mBlkR(4),mBlkR(1):mBlkR(2));
    
    % Define left image keypoint locations using sobel edge detector
    [mX,mY] = meshgrid(mBlkL(1):mBlkL(2),mBlkL(3):mBlkL(4));
    lIn = edge(mI1b);
    mPtsLm = [mX(lIn) mY(lIn)]';
     
    % Make keypoint structure array and extract descriptors for left image
    sKeypoints = cell2struct(num2cell(mPtsLm',2),'pt',2);
    [sKeypoints.size] = deal(0);
    [sKeypoints.angle] = deal(0);
    [sKeypoints.response] = deal(0);
    [sKeypoints.octave] = deal(0);
    [sKeypoints.class_id] = deal(-1);
    mDescriptorsL = descriptors(mI1,sKeypoints);    
    
    % Define right image keypoint locations using sobel edge detector
    [mX,mY] = meshgrid(mBlkR(1):mBlkR(2),mBlkR(3):mBlkR(4));
    lIn = edge(mI2b);
    mPtsRm = [mX(lIn) mY(lIn)]';
     
    % Make keypoint structure array and extract descriptors for right image
    sKeypoints = cell2struct(num2cell(mPtsRm',2),'pt',2);
    [sKeypoints.size] = deal(0);
    [sKeypoints.angle] = deal(0);
    [sKeypoints.response] = deal(0);
    [sKeypoints.octave] = deal(0);
    [sKeypoints.class_id] = deal(-1);
    mDescriptorsR = descriptors(mI2,sKeypoints);

    % Compute matches using hamming distance metric
    [vIdxL,vDist] = knnsearch(double(mDescriptorsL),double(mDescriptorsR), ...
        'NSMethod','exhaustive','Distance','hamming');    
    mPtsLm = mPtsLm(:,vIdxL);
    
    % Save 10 best matches
    [~,vIdx] = sort(vDist);
    try
        mPtsLm = mPtsLm(:,vIdx(1:10));
        mPtsRm = mPtsRm(:,vIdx(1:10));
    catch
    end

    % Save output
    mPtsL = [mPtsL mPtsLm];
    mPtsR = [mPtsR mPtsRm];
    
end

end
end
