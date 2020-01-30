function mZ = meshDenoise(mZ,sInput)

% Parse
vLon = sInput.lon;
vLat = sInput.lat;
vSz = size(mZ);
nullVal = sInput.null;
dT = sInput.params(1);
iN = sInput.params(2);
iV = sInput.params(3);
dBlock = sInput.blockSize;

% Set null values to NaN
if strcmp(nullVal,'dem')
    mZ(mZ < -500 | mZ > 9000) = NaN;
else
    mZ(mZ == nullVal) = NaN;
end

% Mask for null values
lIn = isfinite(mZ);   
             
% Convert to UTM
[mX,mY] = meshgrid(vLon,vLat);
[mPts,iZ,strH] = ll2utm([mX(:) mY(:)],[],[]);
mX = reshape(mPts(:,1),vSz);
mY = reshape(mPts(:,2),vSz);
clear mPts             

% Buffer size
iB = 5;

% Define block vectors
iNumWinX = ceil(vSz(2)/dBlock);
iNumWinY = ceil(vSz(1)/dBlock);
vX = round(linspace(iB+1,vSz(2)-iB,iNumWinX));
vY = round(linspace(iB+1,vSz(1)-iB,iNumWinY));

% Fix vectors if only single block
if length(vX) == 1
    vX = [iB+1 vX];
end
if length(vY) == 1
    vY = [iB+1 vY];
end

% Loops for each block
for iY = 1:length(vY)-1
    for iX = 1:length(vX)-1
        
        % Get point location index for current block
        [mIdxX,mIdxY] = meshgrid(vX(iX)-iB:vX(iX+1)+iB, ...
                                 vY(iY)-iB:vY(iY+1)+iB);
        vIdx = sub2ind(vSz,mIdxY(:),mIdxX(:));
                       
        % Make point matrix               
        mVerts = [mX(vIdx) mY(vIdx) mZ(vIdx)];
        
        % Remove empty points
        mVerts = mVerts(lIn(vIdx),:);
        
        if size(mVerts,1) > 2
        
            % Call main mesh denoising function
            mVertsD = main(mVerts,dT,iN,iV);
            
            % Restore points to original locations
            mVerts = NaN(length(vIdx),3);
            mVerts(lIn(vIdx),:) = mVertsD;
            
            % Remove buffer
            lM = false(size(mIdxX));
            lM(1+iB:end-iB,1+iB:end-iB) = true;
            vIdx = vIdx(lM);
            
            % Save points in matrices
            mX(vIdx) = mVerts(lM,1);
            mY(vIdx) = mVerts(lM,2);
            mZ(vIdx) = mVerts(lM,3);
        
        end
    end
end

% Convert back to geographic coordinate system
mPts = utm2ll([mX(:) mY(:)],iZ,strH);
mX = reshape(mPts(:,1),vSz);
mY = reshape(mPts(:,2),vSz);
clear mPts   
                 
% Compute interpolant
F = scatteredInterpolant(mX(lIn),mY(lIn),mZ(lIn),'linear','none');

% Sample interpolant
[mXs,mYs] = meshgrid(vLon,vLat);
mZ = F(mXs,mYs); 

% Reset null values
mZ(~lIn) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mVerts = main(mVerts,dT,iN,iV)
    
% Compute delaunay triangulation and face normal vectors
sDT = delaunayTriangulation(mVerts(:,1:2));
mFaces = sDT.ConnectivityList;
mNorms = faceNormal(triangulation(mFaces,mVerts));
clear sDT

% Get neighborhoods around faces. Note that the matrix diagonal is included,
% so each face neighborhood includes the face itself. Hence the average of
% the normal vectors around each face includes the normal vector of the
% face. Also, only up to 15 neighbors per face are supported. Any more than
% this consumes too much memory. This may introduce some artifacts in the
% result, especially near edges.
iL = size(mFaces,1);
iP = size(mVerts,1);
mHood = sparse(repmat((1:iL)',1,3),mFaces,1,iL,iP);
mHood = mHood*mHood';
[vR,vC] = find(mHood);
mHood(sub2ind(size(mHood),vR,vC)) = vC;
mHood = sort(mHood,2,'descend');
mHood = full(mHood(:,1:15));
lMask = mHood == 0;
mHood(lMask) = 1;
clear vR vC
    
% Loop for each update to normals
for i = 1:iN
    
    % Get local normals matrix
    mLocalNorms = permute(mNorms(mHood(:),:),[1 3 2]);
    mLocalNorms = reshape(mLocalNorms,[size(mHood) size(mLocalNorms,3)]);    
    mLocalNorms(repmat(lMask,[1 1 3])) = NaN;
    
    % Compute weights matrix
    mNorms = permute(mNorms,[1 3 2]);
    mNorms = repmat(mNorms,[1 size(mLocalNorms,2) 1]);
    mD = dot(mNorms,mLocalNorms,3);
    mH = (mD-dT).^2;
    mH(mD <= dT | isnan(mD)) = NaN;
    
    % Compute normal updates
    mNorms = nansum(repmat(mH,[1 1 3]).*mLocalNorms,2);
    mNorms = permute(mNorms,[1 3 2]);
    mNorms = mNorms./repmat(sqrt(sum(mNorms.^2,2)),1,3);
    
end
clear mD mH

% Get neighborhoods around vertices. Note: up to 10 neighbors are permitted
iL = size(mFaces,1);
iP = size(mVerts,1);
mHood = sparse(repmat((1:iL)',1,3),mFaces,1,iL,iP);
[vR,vC] = find(mHood);
mHood(sub2ind(size(mHood),vR,vC)) = vR;
mHood = sort(mHood,1,'descend')';
mHood = full(mHood(:,1:10));
lMask = mHood == 0;
clear vR vC

% Count number of faces around each vertex (cardinality)
vCard = sum(mHood ~= 0,2);

% Get local normals matrix
mHood(lMask) = 1;
mLocalNorms = permute(mNorms(mHood(:),:),[1 3 2]);
mLocalNorms = reshape(mLocalNorms,[size(mHood) size(mLocalNorms,3)]);
mLocalNorms(repmat(lMask,[1 1 3])) = NaN;

% Loop for each update to the vertices
for i = 1:iV
    
    % Compute face centers
    mCenters = permute(mVerts(mFaces(:),:),[1 3 2]);
    mCenters = reshape(mCenters,[size(mFaces) size(mCenters,3)]);
    mCenters = permute(mean(mCenters,2),[1 3 2]);
    mCenters = permute(mCenters(mHood(:),:),[1 3 2]);
    mCenters = reshape(mCenters,[size(mHood) size(mCenters,3)]);

    % Make vertices matrix
    mVerts3d = permute(mVerts,[1 3 2]);
    mVerts3d = repmat(mVerts3d,[1 size(mLocalNorms,2) 1]);
    
    % Compute vertex updates
    mP = mLocalNorms.*repmat(dot(mLocalNorms,mCenters-mVerts3d,3),[1 1 3]);
    mVerts = mVerts+1./repmat(vCard,1,3).*permute(nansum(mP,2),[1 3 2]);
            
end
end
end
