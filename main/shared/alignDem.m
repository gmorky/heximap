function mT = alignDem(mPts,mDem,lM,vX,vY)

% Moving points and normals
[mX,mY] = meshgrid(vX,vY);
mDemM = points2grid(mPts,vX,vY,'sparse');
lIn = ~isnan(mDemM);
mPtsM = [mX(lIn) mY(lIn) mDemM(lIn)];
mNm = normals(mX,mY,mDemM,true);
mNm = mNm(lIn(:),:);

% Choose subset of moving points
[mPtsM,vIdx] = geoSamplePoints(mPtsM,10000);
mNm = mNm(vIdx,:);

% Remove any remaining NaNs
lIn = all(~isnan(mNm),2);
mNm = mNm(lIn,:);
mPtsM = mPtsM(lIn,:);

% Reference points and normals
lIn = ~isnan(mDem) & ~lM;
mPtsR = [mX(lIn) mY(lIn) mDem(lIn)];
mNr = normals(mX,mY,mDem,true);
mNr = mNr(lIn(:),:);

% Remove any remaining NaNs
lIn = all(~isnan(mNr),2);
mNr = mNr(lIn,:);
mPtsR = mPtsR(lIn,:);

% Parameters
dTol = 0.05;
iMaxIter = 100;
iRejectScale = 3;
iNumLevels = 8;
iSamplingType = 0;
lVis = false;

% Iterative closest point algorithm
mT = icp_mod_point_plane_pyr(mPtsM,mNm,mPtsR,mNr, ...
    dTol,iMaxIter,iRejectScale,[],iNumLevels,iSamplingType,lVis);
