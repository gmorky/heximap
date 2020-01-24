function [mT,sOutput] = geoOptiTrans(mPtsT,strRef,strShpPath,iZ,strH,sOpt)
% Nonlinear optimization of point cloud to match the reference DEM

% Get reference DEM
mBnd = utm2ll([min(mPtsT(1:2,:),[],2) max(mPtsT(1:2,:),[],2)]',iZ,strH);
[mDem,vLon,vLat] = geoGetRefDem(mBnd,strRef);

% Convert reference DEM to points in UTM coordinate system
[mLon,mLat] = meshgrid(vLon,vLat);
mPtsR(:,1:2) = ll2utm([mLon(:) mLat(:)],iZ,strH); clear mLon mLat
mPtsR(:,3) = mDem(:); clear mDem
mPtsR = mPtsR';

% Make reference DEM grid in UTM coordinate system
dX = 90; dY = 90;
vX = floor(min(mPtsR(1,:)))-dX:dX:ceil(max(mPtsR(1,:))+dX);
vY = fliplr(floor(min(mPtsR(2,:)))-dY:dY:ceil(max(mPtsR(2,:))+dY));
mDem = points2grid(mPtsR,vX,vY,'sparse');

% Make unstable terrain grid in UTM coordinate system
lM = polygons2grid(strShpPath,vLon,vLat);
lM = points2grid([mPtsR(1:2,:);double(lM(:)')],vX,vY,'sparse') > 0;
lM = imdilate(lM,strel('disk',3));
clear mPtsR

% Align triangulated points to match reference DEM using ICP algorithm
mT = alignDem(mPtsT,mDem,lM,vX,vY);
mPtsT = [mT * mPtsT; ones(1,size(mPtsT,2))];

% Mask outliers during optimization
lIn = abs(mDem - points2grid(mPtsT,vX,vY,'sparse'));
[mX,mY] = meshgrid(vX,vY);
lIn = interp2(mX,mY,lIn,mPtsT(1,:),mPtsT(2,:)) < 50;

% Optimize orientation of triangulated points to match reference DEM
sOutput = optimizeDem(mPtsT(:,lIn),mDem,lM,vX,vY,sOpt);
