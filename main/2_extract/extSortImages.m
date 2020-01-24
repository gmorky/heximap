function [cHexFile,cM] = extSortImages(cHexFile,cM,hW)
% Determine order in which the Hexagon images were photographed.  This is
% used later when computing stereo disparity maps, to determine which
% images should be considered left or right in each pair.

% Update waitbar
try
waitbar(2/8,hW,'sorting images...')
catch
end

% Get height and width of images
[cH,cW] = cellfun(@(x) size(x,'Image'),cM,'UniformOutput',0);

% Get spatial transformation structures
cT = cellfun(@(x) x.SpatialTrans,cM,'UniformOutput',0);

% Loop for each image
mCenters = zeros(numel(cM),2);
mDir = zeros(numel(cM),2);
for i = 1:numel(cM)
    
    % Define corners. Note y-axis sign is reversed
    mPtsImg = [1 1; cW{i} cH{i}; cW{i} 1; 1 cH{i}];
    mPtsImg(:,2) = -mPtsImg(:,2);
    
    % Transform corners to world coordinates
    mPtsWld = transformPointsForward(cT{i},mPtsImg);

    % Compute image center
    mCenters(i,:) = mean(mPtsWld);
    
    % Compute vector from corner 1 to corner 3
    mDir(i,:) = mPtsWld(3,:) - mPtsWld(1,:);    
    
end

% Determine image order by rotating the images so image x-axes are
% horizontal in the world coordinate system
vPosX = mean([mDir zeros(size(mDir,1),1)]) / ...
   norm(mean([mDir zeros(size(mDir,1),1)]));
mR = vecRotMat(vPosX,[1 0 0]);
mCentersR = (mR * [mCenters zeros(size(mCenters,1),1)]')';
[~,vOrder] = sortrows(mCentersR,1);

% Sort files in order
cHexFile = cHexFile(vOrder);
cM = cM(vOrder);
