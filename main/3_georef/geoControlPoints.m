function [] = geoControlPoints(objL,strRef,strCP,hW)
% Define control points for georeferencing

% Update waitbar
try
waitbar(4/8,hW,'defining control points for initial georeferencing...')
catch
end

% Get control points
switch strCP   
    case 'automatic'
        try
            [mPtsT,mPtsW] = auto(objL,strRef);
        catch
            warning(['Automatic control point selection failed. ' ...
                     'Manual selection is required...'])
            [mPtsT,mPtsW] = manual(objL);
        end
    case 'manual'
        [mPtsT,mPtsW] = manual(objL);
    otherwise
        error('Invalid control point selection parameter.')
end

% Save points in the .mat file
sGeo.Initial.Source = 'self';
sGeo.Initial.WindowTransform = 'none';
sGeo.Initial.GroundControlPoints.Triangulated = mPtsT;
sGeo.Initial.GroundControlPoints.World = mPtsW;
objL.GeorefInfo = sGeo;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mPtsT,mPtsW] = auto(objL,strRef)
    
% Make window in world coordinate system
mWinWld = worldWindow(objL);

% Get reference DEM grid
mBnd = [min(mWinWld);max(mWinWld)];
vB = diff(mBnd)*0.50;
mBnd = mBnd + [-vB;vB];
[mDemR,vLonR,vLatR] = geoGetRefDem(mBnd,strRef);
mDemR = double(mDemR);

% Choose a few random image points and transform to world coordinates
[mPtsI,vIdx] = geoSamplePoints(objL.ImagePoints,100);
mPtsI(1,:) = mPtsI(1,:) + objL.Window(1,1);
mPtsI(2,:) = -(mPtsI(2,:) + objL.Window(1,2));
sSrcInfo = objL.SourceImageInfo;
mPtsW = transformPointsForward(sSrcInfo.SpatialTrans,mPtsI(1:2,:)')';
mPtsW(3,:) = 0;

% Get corresponding triangulated points
mPtsT = objL.TriangulatedPoints;
mPtsT = mPtsT(:,vIdx);

% Use Horn's absolute orientation method to estimate transformation between
% the triangulated points and the world points
lIn = all(isfinite(mPtsT),1) & all(isfinite(mPtsW),1);
sT = absor(mPtsT(1:3,lIn),mPtsW(1:3,lIn),'doScale',1);
mT = sT.M;

% Transform a larger set of triangulated points to world coordinate system
mPtsT = geoSamplePoints(mT*objL.TriangulatedPoints,1E6);

% Make triangulated points grid
dX = mode(abs(diff(vLonR)));
dY = mode(abs(diff(vLatR)));
vLonT = min(mPtsT(1,:)):dX:max(mPtsT(1,:));
vLatT = max(mPtsT(2,:)):-dY:min(mPtsT(2,:));
mDemT = points2grid(mPtsT,vLonT,vLatT,'sparse');
clear mPtsT

% Interpolate holes
mDemRf = inpaint_nans(mDemR);
mDemTf = inpaint_nans(mDemT);

% Choose keypoints
iWinSz = round(mean(size(mDemR))*1/20);
mPtsR = keypoints(mDemRf,iWinSz,5000);
mPtsT = keypoints(mDemTf,iWinSz,500);

% Initialize
vScale = 0.5:0.1:1.5;
iRad = mean(size(mDemT))*1/20;
vRad = round([iRad/2 iRad iRad*2]);
sInliers = struct();

% Loop through scales
for i = 1:length(vScale)
    
    % Extract descriptors
    [mDescriptorsR,mPtsRs] = descriptors(mDemRf,mPtsR,vRad*vScale(i));
    [mDescriptorsT,mPtsTs] = descriptors(mDemTf,mPtsT,vRad);

    % Find potential matches
    [vIdxR,~] = knnsearch(mDescriptorsR',mDescriptorsT');    
    mPtsRs = mPtsRs(:,vIdxR);
    
    % Parameters for RANSAC algorithm
    sRansac.numIter = 2000;
    sRansac.inlierDist = 25*mean(abs([diff(vLonR) diff(vLatR)]));
    
    % Use RANSAC to estimate inliers
    mPtsRw = [vLonR(mPtsRs(1,:)); vLatR(mPtsRs(2,:)); ...
        zeros(1,size(mPtsRs,2))];
    mPtsTw = [vLonT(mPtsTs(1,:)); vLatT(mPtsTs(2,:)); ...
        zeros(1,size(mPtsTs,2))];
    [~,lIn] = estimateTransformRansac(mPtsRw,mPtsTw,sRansac);
    
    % Save output
    sInliers(i).reference = mPtsRw(:,lIn);
    sInliers(i).triangulated = mPtsTw(:,lIn);
    sInliers(i).count = sum(lIn);
  
end

% Choose scale with highest inlier count
[~,iIdx] = max([sInliers.count]);
mPtsW = sInliers(iIdx).reference;
mPtsT = sInliers(iIdx).triangulated;

% Get vertical coordinates for reference DEM points
[mLon,mLat] = meshgrid(vLonR,vLatR);
mPtsW(3,:) = interp2(mLon,mLat,mDemRf,mPtsW(1,:),mPtsW(2,:));
clear mLon mLat

% Get vertical coordinates for triangulated points
mPtsT = mT \ [mPtsT; ones(1,size(mPtsT,2))];
mPts = objL.TriangulatedPoints;
vIdx = knnsearch(mPts(1:2,:)',mPtsT(1:2,:)');
mPtsT = mPts(1:3,vIdx);

% Transpose
mPtsT = mPtsT';
mPtsW = mPtsW';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mPtsT,mPtsW] = manual(objL)

% Make window in world coordinate system
mWinWld = worldWindow(objL);

% Make world coordinate grids
vSz = size(objL,'Image');
[mPts,mY] = meshgrid(linspace(min(mWinWld(:,1)),max(mWinWld(:,1)),vSz(2)), ...
                     linspace(max(mWinWld(:,2)),min(mWinWld(:,2)),vSz(1)));

% Transform world coordinate grids to image coordinates
sSrcInfo = objL.SourceImageInfo;
sT = sSrcInfo.SpatialTrans;
mPts = [mPts(:) mY(:) ones(numel(mPts),1)] / sT.T;
mPts(:,2) = -mPts(:,2);
mPts(:,3) = [];

% Sample image at transformed point locations
mWin = objL.Window;
[mI,mY] = meshgrid(mWin(1):mWin(2),mWin(3):mWin(4));
mPts = interp2(mI,mY,single(objL.ImageFilt),mPts(:,1),mPts(:,2),'linear');
mI = mat2gray(reshape(mPts,vSz));
clear mY mPts

% Show full image
hF1 = figure;
warning('off','images:initSize:adjustingMag');
imshow(adapthisteq(objL.SourceImage)),colormap(bone(256)),hold on
vSz = fliplr(size(objL,'SourceImage'));
text(vSz(2)/2,1,['The orange box indicates the region you '  ...
                   'selected for initial georeferencing.'], ...
                   'Color',[1 0.5 0], ...
                   'Background','k', ...
                   'VerticalAlignment','top');
warning('on','images:initSize:adjustingMag');
   
% Rectangle for chosen window
dScale = 10;
mWin = objL.Window/dScale;
vRec = [min(mWin(:,1)) min(mWin(:,2)) ...
    abs(diff(mWin(:,1))) abs(diff(mWin(:,2)))];
rectangle('Position',vRec,'EdgeColor',[1 0.5 0])

%Show (roughly georeferenced) image window
hF2 = figure;
warning('off','images:initSize:adjustingMag');
imshow(mI,'XData',[min(mWinWld(:,1)) max(mWinWld(:,1))], ...
          'YData',[max(mWinWld(:,2)) min(mWinWld(:,2))])
set(gca,'YDir','normal')
axis equal,colormap(bone(256)),hold on
warning('on','images:initSize:adjustingMag');
set(gcf,'name','Control points','numbertitle','off')
clear mI

%Initialize variables
iC = 1;
iNumPts = 5;
mPtsI = zeros(iNumPts,2);
mPtsW = zeros(iNumPts,3);

%Loop for each point
while iC <= iNumPts
    
    % Error if user closes figure window
    if ~ishandle(hF2)
        error('Process aborted by user.')
    end
    
    % Update text
    hT = text(mean(mWinWld(:,1)),max(mWinWld(:,2)), ...
        ['Click and drag to place control point ' ...
        num2str(iC) ' of ' num2str(iNumPts) ...
        '. Double click on point when finished.'], ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','Bottom', ...
        'Color',[1 0.5 0], ...
        'Background','k');

    %User inputs image point
    hP = impoint; zoom out
    vPt = wait(hP);
    
    % Initialize
    lValidInput = false(1,3);
    lCancelled = false;
    
    % Make sure user input is valid
    while any(~lValidInput)

        %User inputs world point
        cP = {'Control Point Longitude (dec):', ...
              'Control Point Latitude (dec):', ...
              'Elevation (meters):'};
        strTitle = 'Input';
        mSz = repmat([1 45],3,1);       
        cIn = inputdlg(cP,strTitle,mSz);

        % If user hit cancel button
        if numel(cIn) == 0
            lValidInput = true(1,3);
            lCancelled = true;
            continue
        end
        
        % Convert strings to numeric
        vIn = str2double(cIn);

        % Check validity of longitude
        if abs(vIn(1)) <= 180
            lValidInput(1) = true;
        end
        
        % Check validity of latitude
        if abs(vIn(2)) <= 90
            lValidInput(2) = true;
        end
        
        % Check validity of elevation
        if vIn(3) > -500 && vIn(3) < 9000
            lValidInput(3) = true;
        end
        
        % Continue if all input is valid
        if all(lValidInput)           
            continue
        end

        % Warn
        uiwait(warndlg('    Invalid input. Please try again.','Warning'));
        
        % Reset
        lValidInput = false(1,3);

    end
    
    % Delete text and control point handles
    delete(hT),delete(hP)
    
    if ~lCancelled
    
        % Plot the selected point
        scatter(vPt(1),vPt(2),'s','MarkerEdgeColor',[1 0.5 0], ...
            'LineWidth',1),pause(0.1)
        
        %Save points
        mPtsI(iC,:) = transformPointsInverse(sT,vPt) .* [1 -1];
        mPtsW(iC,:) = str2double(cIn)';
        
        % Increment count
        iC = iC + 1;
        
    end
end

% Close the figure windows
close(hF1),close(hF2),pause(0.1)

% Find triangulated points corresponding to the selected image points
mPts = objL.ImagePoints';
vIdx = dsearchn(mPts(:,1:2),mPtsI);
mPts = objL.TriangulatedPoints';
mPtsT = mPts(vIdx,1:3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mWinWld = worldWindow(objL)
    
% Get data
sSrcInfo = objL.SourceImageInfo;
sT = sSrcInfo.SpatialTrans;
mWin = objL.Window;

% Make window in world coordinate system
vUL = [min(mWin(:,1)) min(mWin(:,2))];
vUR = [max(mWin(:,1)) min(mWin(:,2))];
vLR = [max(mWin(:,1)) max(mWin(:,2))];
vLL = [min(mWin(:,1)) max(mWin(:,2))];
mWin = [vUL;vUR;vLR;vLL];
mWin(:,2) = -mWin(:,2);
mWinWld = transformPointsForward(sT,mWin);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mPtsOut = keypoints(mGrid,iWinSz,iNumPts)
% Choose candidate points by finding local minimum and maximum values of
% image

% Initialize
mIdx = reshape(1:numel(mGrid),size(mGrid));
iInc = round(iWinSz/2);
iEmpty = 0;

% Loop through search window sizes
for iW = 1:round(min(size(mGrid))/iInc)
    
    % Increment search window size and reset loop count
    iWinSz = iWinSz + iInc;
    iCount = 1;
    
    % Loop through search window locations
    for iR = 1:round(iWinSz/4):size(mGrid,1) - iWinSz
        for iC = 1:round(iWinSz/4):size(mGrid,2) - iWinSz
            
            % Extract sub image and sub coordinates
            mSubI = mGrid(iR:iR+iWinSz-1,iC:iC+iWinSz-1);
            mSubIdx = mIdx(iR:iR+iWinSz-1,iC:iC+iWinSz-1);
            
            if sum(mSubI == iEmpty) <= sum(mSubI ~= iEmpty)
                
                % Find minimum values and their locations
                [dVal,iIdx] = min(mSubI(:));
                mMinIdx(iW,iCount) = mSubIdx(iIdx);
                mMinVal(iW,iCount) = dVal;
                
                % Find maximum values and their locations
                [dVal,iIdx] = max(mSubI(:));
                mMaxIdx(iW,iCount) = mSubIdx(iIdx);
                mMaxVal(iW,iCount) = dVal;
                
                % Increment loop count
                iCount = iCount + 1;
                
            end
        end
    end
end

% Count number of duplicate elements in index matrices
vMaxIdxU = unique(mMaxIdx(mMaxVal ~= iEmpty));
vMinIdxU = unique(mMinIdx(mMinVal ~= iEmpty));
vMaxCnt = hist(mMaxIdx(:),vMaxIdxU)';
vMinCnt = hist(mMinIdx(:),vMinIdxU)';

% Sort the matrices, keep the points with the most counts, convert to
% rows and columns
mMax = flipud(sortrows([vMaxIdxU vMaxCnt],2));
mMin = flipud(sortrows([vMinIdxU vMinCnt],2));
try
mMax = mMax(1:iNumPts,:);
mMin = mMin(1:iNumPts,:);
catch
end
[vYmax,vXmax] = ind2sub(size(mGrid),mMax(:,1));
[vYmin,vXmin] = ind2sub(size(mGrid),mMin(:,1));

% Assign output
mPtsMax = [vXmax vYmax];
mPtsMin = [vXmin vYmin];
mPtsOut = [mPtsMax;mPtsMin]';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mDescriptors,mPts] = descriptors(mGrid,mPts,vRad)
    
% Initialize
mDescriptors = [];

% Loop through radii
for i = 1:length(vRad)

    % Parameterize the circle
    vT = linspace(0,2*pi,360);
    mH = repmat(mPts(1,:),[length(vT) 1]);
    mK = repmat(mPts(2,:),[length(vT) 1]);
    mX = mH + vRad(i)*cos(repmat(vT(:),[1 size(mPts,2)]));
    mY = mK + vRad(i)*sin(repmat(vT(:),[1 size(mPts,2)]));

    % Sample the grid
    mZ = interp2(double(mGrid),mX,mY);

    % Shift columns so max pixel values are first
    [iR,iC] = size(mZ);
    [~,vShiftIdx] = max(mZ);
    mS = full(sparse(mod(-vShiftIdx+1,iR)+1,1:iC,1,iR,iC));
    mZ = ifft(fft(mZ).*fft(mS),'symmetric');

    % Divide by standard deviation to normalize, then subtract mean to
    % center
    mZ = mZ ./ repmat(std(mZ),[size(mZ,1) 1]);
    mZ = mZ - repmat(mean(mZ),[size(mZ,1) 1]);
    
    % Save output
    mDescriptors = [mDescriptors;mZ];

end

% Remove points near edges
lIn = all(~isnan(mDescriptors),1);
mDescriptors = mDescriptors(:,lIn);
mPts = mPts(:,lIn);

end
end