function [] = geoOptimize(cL,iWinIdx,strRef,strShpPath,lVis,lPoly,hW)
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
    % least 1/3 of the image points were successfully triangulated. If not,
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
    while dRMSE > 20 && iCount <= iAttempts
        
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
        
        % Apply transformations from alignment and optimization routines
        iIdx = vWinIdx(iW);
        mPtsT = cL{iIdx}.TriangulatedPointsGeoref;
        for iT = 1:numel(cOutput)
            mPtsT = [cT{iT} * mPtsT; ones(1,size(mPtsT,2))];
            mPtsT = transformUsingSolverVar(mPtsT,cOutput{iT});
        end

        % Convert from UTM to WGS84 latitude and longitude
        mPtsT = utm2ll(mPtsT',iZ,strH)';
        
        % Get boundaries then save points
        mBnd = [min(mPtsT(1:2,:),[],2) max(mPtsT(1:2,:),[],2)]';
        cL{iIdx}.TriangulatedPointsGeoref = mPtsT;
        clear mPtsT
        
        % Save reference DEM     
        [cL{iIdx}.ReferenceDem,cL{iIdx}.ReferenceLon, ...
            cL{iIdx}.ReferenceLat] = geoGetRefDem(mBnd,strRef);

        % Save georeferencing transformation data
        sGeo = cL{iIdx}.GeorefInfo;
        sGeo.Final.AlignmentOutput = cT;
        sGeo.Final.OptimizationOutput = cOutput;
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
