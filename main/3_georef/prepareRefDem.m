function [vX,vY,mDem,lM] = prepareRefDem(strRef,strShpPath,mBnd,iZ,strH)

% Get reference DEM
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
