function [mT,lKeep] = estimateTransformRansac(mPts1,mPts2,sRansac)
% Estimate 2D transformation using simple RANSAC implementation

% Check for errors
if isempty(mPts1) || isempty(mPts2) || any(size(mPts1) ~= size(mPts2))
    error('Point arrays must be non-empty and the same size.')
end
if size(mPts1,1) ~= 3 || size(mPts2,1) ~= 3
    error('Point arrays must be 3xN in size.')
end

% Compute transformation using ransac
[mT,lKeep] = ransac(mPts1,mPts2,sRansac);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mT,lIn] = ransac(mPts1,mPts2,sRansac)

% Check if user specified consensus metric function handle
lDefault = true;
if isfield(sRansac,'consensusMetric')
    if isa(sRansac.consensusMetric,'function_handle')
        lDefault = false;
        hFun = sRansac.consensusMetric;
    end
end

% Loop for each iteration
iSetSize = 3;
vC = zeros(1,sRansac.numIter);
cT = cell(1,sRansac.numIter);
warning('off','images:maketform:conditionNumberofAIsHigh')
warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:singularMatrix')
for i = 1:sRansac.numIter

    % Fit transformation to random subset of data points
    vIdx = randperm(size(mPts1,2),min([iSetSize size(mPts1,2)]));
    mT = trans(mPts1(:,vIdx),mPts2(:,vIdx));
    
    if lDefault
        
        % Count number of data points in the consensus set
        vC(i) = sum(transDist(mT,mPts1,mPts2) < sRansac.inlierDist);
        
    else
        
        % Call user-defined consensus metric function
        sInput.trans = mT;
        sInput.points1 = mPts1;
        sInput.points2 = mPts2;
        vC(i) = hFun(sInput);
        
    end
    
    % Save transformation
    cT{i} = mT;
   
end
warning('on','images:maketform:conditionNumberofAIsHigh')
warning('on','MATLAB:nearlySingularMatrix')
warning('on','MATLAB:singularMatrix')

% Choose transformation with the most inliers
[~,iIdx] = max(vC);
mT = cT{iIdx};
lIn = transDist(mT,mPts1,mPts2) < sRansac.inlierDist;

% Estimate final transformation using entire consensus set
mT = trans(mPts1(:,lIn),mPts2(:,lIn));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mT = trans(mPts1,mPts2)
    
% Use Horn's absolute orientation method to compute initial transformation
lIn = all(isfinite(mPts1),1) & all(isfinite(mPts2),1);
sT = absor(mPts2(:,lIn),mPts1(:,lIn),'doScale',1);
mT = sT.M;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vDist = transDist(mT,mPts1,mPts2)    
    
% Transform the points
mPts2 = mT * [mPts2; ones(1,size(mPts2,2))];

% Compute distance between point pairs
vDist = sqrt(sum((mPts2(1:3,:) - mPts1(1:3,:)).^2));

end
end
