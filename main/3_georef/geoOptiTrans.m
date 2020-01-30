function [mT,sOutput] = geoOptiTrans(mPtsT,strRef,strShpPath,iZ,strH,sOpt)
% Nonlinear optimization of point cloud to match the reference DEM

% Boundaries
mBnd = utm2ll([min(mPtsT(1:2,:),[],2) max(mPtsT(1:2,:),[],2)]',iZ,strH);

% Prepare reference DEM
[vX,vY,mDem,lM] = prepareRefDem(strRef,strShpPath,mBnd,iZ,strH);

% Align triangulated points to match reference DEM using ICP algorithm
mT = alignDem(mPtsT,mDem,lM,vX,vY);
mPtsT = [mT * mPtsT; ones(1,size(mPtsT,2))];

% Mask outliers during optimization
[mX,mY] = meshgrid(vX,vY);
lIn = abs(interp2(mX,mY,mDem-points2grid(mPtsT,vX,vY,'sparse'), ...
    mPtsT(1,:),mPtsT(2,:))) < 50;

% Optimize orientation of triangulated points to match reference DEM
sOutput = optimizeDem(mPtsT(:,lIn),mDem,lM,vX,vY,sOpt);
