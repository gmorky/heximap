function sOutput = optimizeDem(mPts,mDem,lM,vX,vY,sOpt)

% Initialize variables for solver
vVar = [zeros(1,6) ones(1,4)];
iMaxIter = sOpt.maxIterations;
vC = mean(mPts(1:3,:),2);
lVis = sOpt.visualize;
lPoly = sOpt.polySurf;

% Global variables for saving across functions
dError = NaN; dPrevError = Inf; 
vVarCopy = []; vFinalVar = [];
sPoly = 'none'; sFinalPoly = 'none';
dRMSE = NaN; dFinalRMSE = NaN; vRMSE = NaN(1,iMaxIter);

% Create matrix of upper and lower boundaries for scaling
vBnd = [sOpt.rotation sOpt.translation sOpt.scale sOpt.globalScale];
mBnd = repmat(vVar,2,1) + [-vBnd;vBnd];

% Remove variables which won't be included during optimization.  Index
% vector specifies variable locations in array
lIn = vBnd > 0;
vVar = vVar(lIn);
mBnd = mBnd(:,lIn);
vIdx = 1:length(vBnd);
vIdx = vIdx(lIn);

% Scale the variables
for i = 1:length(vVar)
    vVar(i) = (vVar(i) - mBnd(1,i)) / (mBnd(2,i) - mBnd(1,i));
end

% Open figure window if user wants to see animation
if lVis == 1
    hF = figure;
end

% Minimize the objective function using downhill simplex solver
sOpt = optimset('Display','off','MaxIter',iMaxIter,'OutputFcn',@iter);
fminsearch(@(vVar) optError(mPts,mDem,lM,vX,vY,mBnd,vC,vIdx, ...
    lVis,lPoly,vVar),vVar,sOpt);

% Close figure window if necessary
if lVis == 1
    close(hF)
end

% Output structure
sOutput.scaleBounds = mBnd;
sOutput.rotationCenter = vC;
sOutput.variablesIndex = vIdx;
sOutput.variables = vFinalVar;
sOutput.polySurf = sFinalPoly;
sOutput.verticalRMSE = dFinalRMSE;
sOutput.direction = 'forward';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dOutput = optError(mPts,mDemR,lM,vX,vY,mBnd,vC,vIdx,lVis, ...
        lPoly,vVar)
% Error function for downhill simplex solver

% Save copy of current variables
vVarCopy = vVar;

% Apply transformation
sInput.scaleBounds = mBnd;
sInput.rotationCenter = vC;
sInput.variablesIndex = vIdx;
sInput.variables = vVar;
sInput.direction = 'forward';
mPts = transformUsingSolverVar(mPts,sInput);

% Make grid for moving points
mDemT = points2grid(mPts,vX,vY,'sparse');

% Take difference of DEMs and mask unstable terrain
mD = mDemR - mDemT;
mD(lM) = NaN;

% Remove outliers
vQ = quantile(mD(:),[0.05 0.95]);
mD(mD < vQ(1) | mD > vQ(2)) = NaN;

if lPoly

    % Fit 3rd order polynomial surface to residual elevations
    lValid = ~isnan(mD);
    [mX,mY] = meshgrid(vX,vY);
    strTerms = 'x,y,x*y,x^2,y^2,x^3,x*x*y,x*y*y,y^3';
    sPoly = polyfitn([mX(lValid) mY(lValid)],mD(lValid),strTerms);
    mSurf = reshape(polyvaln(sPoly,[mX(:) mY(:)]),size(mD));

    % Apply the correction
    mD = mD - mSurf;
    mDemT = mDemT + mSurf;

end

% Compute error
dError = sum(mD(~isnan(mD)).^2);

% Compute vertical root-mean-square error
dRMSE = sqrt(sum((mD(~isnan(mD))).^2)/numel(mD(~isnan(mD))));

% Visualization
if lVis
    try
    visualize(mDemR,lM,mDemT,mD);
    catch
    end
end

% Objective function output
dOutput = dError;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lStop = iter(vVar,sValues,strState)

% Iteration number
iIter = sValues.iteration;
    
% Initialize stop flag
lStop = false;

% Save current RMSE value in vector
if iIter > 0
    vRMSE(iIter) = dRMSE;
end

% If the previous 10 iterations all resulted in similar RMSE values, stop
% the solver
if iIter > 10
    vTest = vRMSE(iIter-10:iIter-1);
    if (max(vTest) - min(vTest)) < 0.1
        lStop = true;
    end
end

% Save output parameters        
if dError < dPrevError     
        dPrevError = dError;
        vFinalVar = vVarCopy;
        sFinalPoly = sPoly;
        dFinalRMSE = dRMSE;
end

% Display vertical root-mean-square error   
if strcmpi(strState,'iter')       
    disp(['  iteration ' num2str(iIter) ...
        '... approximate relative vertical RMSE (m): ' num2str(dRMSE)])  
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = visualize(mDemR,lM,mDemT,mD)
        
% Compute histogram of pixel differences
[vC2,dN2] = hist(mD(~lM & ~isnan(mD)),256);   

% Show the images
set(0,'currentfigure',hF);
clf,colordef(gcf,'none')
vR = [min(mDemR(:)) max(mDemR(:))];
subplot(2,2,1), hI = imagesc(mDemT,vR);
mAlpha = zeros(size(mDemT));
mAlpha(~isnan(mDemT)) = 1;
set(hI,'AlphaData',mAlpha)
colormap(colormapDEM)
axis off
title('Moving DEM')
subplot(2,2,2), hI = imagesc(mDemR) ;
mAlpha = zeros(size(mDemR));
mAlpha(~isnan(mDemR)) = 1;
set(hI, 'AlphaData', mAlpha)
axis off
title('Reference DEM')
subplot(2,2,3), hI = imagesc(mD);
mAlpha = zeros(size(mD));
mAlpha(~isnan(mD)) = 1;
set(hI, 'AlphaData',mAlpha)
axis off
title('Difference')
subplot(2,2,4), bar(dN2,vC2)
xlabel('Elevation Difference (m)'),ylabel('Count')
title('Histogram')
xlim([-100 100])
drawnow

end
end