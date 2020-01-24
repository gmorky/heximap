function vShift = shiftDem(mPts,mDem,lM,vX,vY,vScale)

% Make grid for moving points
mDemT = points2grid(mPts,vX,vY,'sparse'); clear mPts

% Compute slope and aspect of reference DEM
[mP,mQ] = gradient(mDem,vX*vScale(1),vY*vScale(2));
mS = atand(sqrt(mP.^2+mQ.^2));
mA = atan2d(mP,mQ) + 180;
clear mGx mGy

% Compute elevation difference
mD = mDem - mDemT;

% Exclude unstable terrain
mD(lM) = NaN;

% Exclude outlier elevations
vQ = quantile(mD(:),[0.02 0.98]);
mD(mD < vQ(1) | mD > vQ(2)) = NaN;

% Exclude extreme slopes and missing pixels
lIn = mS > 5 & mS < 75 & ~isnan(mD);

% Fit cosine equation to the data
F = @(x,xdata) x(1)*cosd(x(2)-xdata)+x(3);
sOpt = optimset('Display','off');
vOutput = lsqcurvefit(F,[0 0 0],mA(lIn),mD(lIn)./tand(mS(lIn)),[],[],sOpt);

% Compute shift vector
vShift(1) = vOutput(1)*sind(vOutput(2));
vShift(2) = vOutput(1)*cosd(vOutput(2));
vShift(3) = vOutput(3)*tand(mean(mS(lIn)));

% Divide by scale factor
vShift(1:2) = vShift(1:2) ./ vScale;
