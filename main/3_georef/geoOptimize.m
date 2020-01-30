function [] = geoOptimize(cL,iWinIdx,strRef,strShpPath,lVis,lPoly, ...
    lWin,dMaxDelZ,hW)
% Nonlinear optimization of all windows to match the reference DEM

% Get UTM zone and hemisphere
sGeo = cL{iWinIdx}.GeorefInfo;
sT = sGeo.Initial.Triangulated2WorldTransform;
iZ = sT.zone;
strH = sT.hemi;

% Get regions
vROI = cell2mat(cellfun(@(x) x.RegionID,cL,'Uni',0))';
vUniqueROI = unique(vROI);
iNumROI = length(vUniqueROI);
       
% Parameters for nonlinear optimization
sOpt.rotation = [5 5 5];
sOpt.translation = [5E3 5E3 5E3];
sOpt.scale = [1 1 1];
sOpt.globalScale = 2;
sOpt.maxIterations = 200;
sOpt.visualize = lVis;
sOpt.polySurf = lPoly;

% Loop through each region
for iR = 1:iNumROI
    
    try
    
    % Update waitbar and command window
    try
    waitbar(7/8,hW,['region ' num2str(iR) ' of ' num2str(iNumROI) ':' ...
        ' optimizing point cloud orientation...'])
    catch
    end   
    disp(['optimizing region ' num2str(iR) ' of ' num2str(iNumROI) '...'])
    
    % Find windows belonging to current ROI
    vWinIdx = find(vROI == vUniqueROI(iR));

    % Loop through each window in current region, and determine whether at
    % least 33% of the image points were successfully triangulated. If not,
    % don't use points from the window during optimization.
    lSample = true(size(vWinIdx));
    iCount = 1;
    for iW = 1:numel(vWinIdx)
        iIdx = vWinIdx(iW);
        [~,iPointCount] = size(cL{iIdx},'TriangulatedPoints');
        if iPointCount / ...
                (prod(size(cL{iIdx},'Image')./cL{iIdx}.DisparityScale)) ...
                < 0.33
            lSample(iCount) = false;
        end 
        iCount = iCount + 1;
    end
    
    % Initialize
    iNumPts = round(1E6/sum(lSample));
    mPtsT = [];
    
    % Loop through each window in current region
    vWinIdxS = vWinIdx(lSample);
    for iW = 1:numel(vWinIdxS)
    
        % Choose random subset of triangulated points (to conserve memory)
        iIdx = vWinIdxS(iW);
        mPts = geoSamplePoints(cL{iIdx}.TriangulatedPointsGeoref,iNumPts);
        mPtsT = [mPtsT mPts];
        clear mPts
    
    end
    
    % Nonlinear optimization
    iCount = 1; dRMSE = Inf; iAttempts = 3; cT = {}; cOutput = {};
    while dRMSE > 25 && iCount <= iAttempts
        
        % Update command window
        disp(['  solver attempt ' num2str(iCount) ': '])
        
        % Optimize
        [cT{iCount},cOutput{iCount}] = geoOptiTrans( ...
            mPtsT,strRef,strShpPath,iZ,strH,sOpt);
        
        % Apply transformations
        mPtsT = [cT{iCount} * mPtsT; ones(1,size(mPtsT,2))];
        mPtsT = transformUsingSolverVar(mPtsT,cOutput{iCount});     
        
        % Save RMSE and increment count
        dRMSE = cOutput{iCount}.verticalRMSE;
        iCount = iCount + 1;
        
    end
    clear mPtsT
    
    % Loop through each window in current region
    for iW = 1:numel(vWinIdx)
        
        % Apply transformations from previous alignment and optimization
        % routines
        iIdx = vWinIdx(iW);
        mPtsT = cL{iIdx}.TriangulatedPointsGeoref;
        for iT = 1:numel(cOutput)
            mPtsT = [cT{iT} * mPtsT; ones(1,size(mPtsT,2))];
            mPtsT = transformUsingSolverVar(mPtsT,cOutput{iT});
        end
        
        % Update command window
        disp(['optimizing window ' num2str(iW) ' of region ' ...
            num2str(iR) '...'])
        
        % Additional optimization of window if data coverage is at least
        % 33%
        mBnd = utm2ll([min(mPtsT(1:2,:),[],2) max(mPtsT(1:2,:),[],2)]', ...
            iZ,strH);
        [vX,vY,~,lM] = prepareRefDem(strRef,strShpPath,mBnd,iZ,strH);
        [mX,mY] = meshgrid(vX,vY);
        lIn = true(1,size(mPtsT,2));
        lIn(interp2(mX,mY,double(lM),mPtsT(1,:),mPtsT(2,:)) > 0) = false;
        lW = false;
        if (sum(lIn) / ...
                (prod(size(cL{iIdx},'Image')./cL{iIdx}.DisparityScale)) ...
                > 0.33) && lWin
            [cT{end+1},cOutput{end+1}] = ...
                geoOptiTrans(geoSamplePoints(mPtsT,iNumPts),strRef, ...
                strShpPath,iZ,strH,sOpt);
            mPtsT = [cT{end} * mPtsT; ones(1,size(mPtsT,2))];
            mPtsT = transformUsingSolverVar(mPtsT,cOutput{end});
            lW = true;
        end

        % Convert from UTM to WGS84 latitude and longitude
        mPtsT = utm2ll(mPtsT',iZ,strH)';
        
        % Remove any points which exceed the elevation difference threshold
        mBnd = [min(mPtsT(1:2,:),[],2) max(mPtsT(1:2,:),[],2)]';
        [mDem,vLon,vLat] = geoGetRefDem(mBnd,strRef);
        [mLon,mLat] = meshgrid(vLon,vLat);
        lCut = abs(mPtsT(3,:) - ...
            interp2(mLon,mLat,mDem,mPtsT(1,:),mPtsT(2,:))) > dMaxDelZ;
        clear mLon mLat
        mPtsT(:,lCut) = NaN;
        
        % Save triangulated points
        cL{iIdx}.TriangulatedPointsGeoref = mPtsT;
        clear mPtsT
        
        % Save reference DEM
        cL{iIdx}.ReferenceDem = mDem;
        cL{iIdx}.ReferenceLon = vLon;
        cL{iIdx}.ReferenceLat = vLat;

        % Save georeferencing transformation data
        sGeo = cL{iIdx}.GeorefInfo;
        sGeo.Final.AlignmentOutput = cT;
        sGeo.Final.OptimizationOutput = cOutput;
        sGeo.Final.WindowOptimizedSeparately = lW;
        cL{iIdx}.GeorefInfo = sGeo;

        % Save vertical root mean square error
        sAccuracy = cL{iIdx}.Accuracy;
        sAccuracy.VerticalRMSError = cOutput{end}.verticalRMSE;
        cL{iIdx}.Accuracy = sAccuracy;
        
    end
    
    catch objExc

    warning(objExc.message)
    warning(['Optimization of region ' num2str(iR) ' failed, skipping...'])
    
    end   
end
