function [] = extBundleAdjust(strHexPath,cL,cR,hW,cReg)
% Minimize reprojection error using nonlinear least squares solver

% Update waitbar
try
set(get(findobj(hW,'type','axes'),'title'), 'string', ...
    ['region ' cReg{1} ' of ' cReg{2} ': performing bundle adjustment...'])
pause(0.1)
catch
end

% Choose matrices for initial solver guess
try
    [mK1,mK2,mP1,mP2] = chooseInitialMatrices(strHexPath);
catch
    mK1 = cL{1}.IntrinsicMatrix;
    mK2 = cR{1}.IntrinsicMatrix;
    mP1 = cL{1}.PoseMatrix;
    mP2 = cR{1}.PoseMatrix;
end

% Initialize
dFval = Inf; iCount = 1;

% Perform bundle adjustment 3 times, using output from previous solver run
% as initial guess. This helps if convergence isn't reached in
% single run
while dFval > 100 && iCount <= 3
    
    % Perform bundle adjustment
    [mK1,mK2,mP1,mP2,vError,dFval] = bundleAdjust(cL,cR,mK1,mK2,mP1,mP2);

    % Increment count
    iCount = iCount + 1;
    
end

% Save output
saveOutput(cL,cR,mK1,mK2,mP1,mP2,vError,dFval);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [mK1,mK2,mP1,mP2] = chooseInitialMatrices(strHexPath)
% Choose pose and intrinsic matrices with highest accuracy for initial
% solver guess (if previous matrices have been computed)

% Initialize
cFileL = getFiles(strHexPath,'Left.mat');
cFileR = getFiles(strHexPath,'Right.mat');
dAccuracy = Inf;

% Loop through mat files
for iW = 1:numel(cFileL)
    try
        
        % Get previously computed bundle adjust accuracy
        objL = matfile(cFileL{iW});
        objR = matfile(cFileR{iW});
        sA = objL.Accuracy;
        
        % Test accuracy against previous low
        if sA.BundleAdjust < dAccuracy
            
            % Save matrices
            mP1 = objL.PoseMatrix;
            mP2 = objR.PoseMatrix;
            mK1 = objL.IntrinsicMatrix;
            mK2 = objR.IntrinsicMatrix;
            
            % Save current accuracy value as new test
            dAccuracy = sA.BundleAdjust;
            
        end
        
    catch
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mK1,mK2,mP1,mP2,vError,dFval] = bundleAdjust(cL,cR,mK1,mK2, ...
        mP1,mP2)

% Choose 2000 random points from the image windows
iNumWin = numel(cL);
iNumPts = round(2000/iNumWin);
mPts1 = []; mPts2 = [];
for i = 1:iNumWin
    
    % Get points
    mPts = cL{i}.PointMatches;
    mHexWinL = cL{i}.Window;
    mPts(1,:) = mPts(1,:) + mHexWinL(1) - 1;
    mPts(2,:) = mPts(2,:) + mHexWinL(3) - 1;
    vIdx = randperm(length(mPts),min([iNumPts length(mPts)]));
    mPts1 = [mPts1 mPts(:,vIdx)];
    mPts = cR{i}.PointMatches;
    mHexWinR = cR{i}.Window;
    mPts(1,:) = mPts(1,:) + mHexWinR(1) - 1;
    mPts(2,:) = mPts(2,:) + mHexWinR(3) - 1;
    mPts2 = [mPts2 mPts(:,vIdx)];
    
end
mPts1 = [mPts1;ones(1,length(mPts1))];
mPts2 = [mPts2;ones(1,length(mPts2))];
clear mPts

% Extract focal length, principle point, camera pose angles and positions
vF = [mK1(1,1);mK1(2,2)];
vPP1 = mK1(1:2,3);
vPP2 = mK2(1:2,3);
vPose2 = [flipud(rotro2eu('zyx', mP2(1:3,1:3))) * 180/pi; mP2(1:3,4)];

% Assign variables to vector
vVar = [vF;vPP1;vPP2;vPose2];

% Define scaling
dFx = 5; dFy = 5;
dPPx = 5; dPPy = 5;
dRx = 5; dRy = 5; dRz = 5;
dTx = abs(mP1(1,4)-mP2(1,4))/2; 
dTy = abs(mP1(2,4)-mP2(2,4))/2; 
dTz = abs(mP1(3,4)-mP2(3,4))/2;

% Create matrix of upper and lower boundaries
vBnd = [dFx dFy repmat([dPPx dPPy],1,2) dRx dRy dRz dTx dTy dTz]';
mBnd = repmat(vVar,1,2) + [-vBnd vBnd];

% Scale the parameters
for i = 1:length(vVar)
    vVar(i) = (vVar(i) - mBnd(i,1)) / (mBnd(i,2) - mBnd(i,1));
end

% Specify lower and upper boundaries for solver
vLB = zeros(size(vVar));
vUB =  ones(size(vVar));

% Minimize reprojection error using non-linear least squares solver
sOpt = optimset('Display', 'iter', 'MaxIter', 100, 'TolX', 1E-6, ...
    'TolFun', 1E-6);
[vVar,dFval] = lsqnonlin(@(vVar) reprojectionError(mPts1,mPts2, ...
    mP1,mBnd,vVar),vVar,vLB,vUB,sOpt);

% Compute optimized intrinsic and pose matrices
[mK1,mK2,mP2] = computeMatrices(mBnd,vVar);

% Compute final reprojection error
[~,vError] = triangulate(mPts1,mPts2,mP1,mP2,mK1,mK2,true);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vError = reprojectionError(mPts1,mPts2,mP1,mBnd,vVar)
% Compute reprojection error by triangulating the points in 3D space, then
% projecting back onto the image planes.

% Compute intrinsic and pose matrices
[mK1o,mK2o,mP2o] = computeMatrices(mBnd,vVar);

% Triangulate 3D points and output the reprojection error
[~,vError] = triangulate(mPts1,mPts2,mP1,mP2o,mK1o,mK2o,true);

% Remove outliers
vError(vError >= quantile(vError,0.98)) = 0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mK1,mK2,mP2] = computeMatrices(mBnd,vVar)
% Compute intrinsic and pose matrices

% Descale the parameters
for iR = 1:length(vVar)
    vVar(iR) = mBnd(iR,1) + vVar(iR) * (mBnd(iR,2) - mBnd(iR,1));
end

% Define the first camera intrinsic matrix
mK1 = [vVar(1) 0          vVar(3)
       0          vVar(2) vVar(4)
       0          0          1];

% Define the second camera intrinsic matrix
mK2 = [vVar(1) 0          vVar(5)
       0          vVar(2) vVar(6)
       0          0          1];

% Make pose matrix
mR2 = makeRotMat(vVar(7),vVar(8),vVar(9),'deg');
vT2 = [vVar(10) vVar(11) vVar(12)]';
mP2 = [mR2 vT2];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = saveOutput(cL,cR,mK1,mK2,mP1,mP2,vError,dFval)
        
% Loop through each window
for iW = 1:numel(cL)

    % Save output
    cL{iW}.IntrinsicMatrix = mK1;
    cR{iW}.IntrinsicMatrix = mK2;
    cL{iW}.PoseMatrix = mP1;
    cR{iW}.PoseMatrix = mP2;

    % Save accuracy stats
    sAccuracy = cL{iW}.Accuracy;
    sAccuracy.Reprojection = vError;
    sAccuracy.BundleAdjust = dFval;
    cL{iW}.Accuracy = sAccuracy;

end

end
end