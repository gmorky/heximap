function [] = geoApplyTrans(cL,cR,iWinIdx,hW)
% Apply initial georeferencing transformations

% Get transformations from chosen window
sGeo = cL{iWinIdx}.GeorefInfo;
sT = sGeo.Initial.Triangulated2WorldTransform;
mTa = sGeo.Initial.AlignmentOutput;
sOutput = sGeo.Initial.OptimizationOutput;

% Get regions
vROI = cell2mat(cellfun(@(x) x.RegionID,cL,'Uni',0))';
vUniqueROI = unique(vROI);
iNumROI = length(vUniqueROI);

% Loop through each region
for iR = 1:iNumROI
    
    try
    
    % Update waitbar and command window
    try
    waitbar(6/8,hW,['region ' num2str(iR) ' of ' num2str(iNumROI) ':' ...
        ' applying initial transformations...'])
    catch
    end   
    disp(['transforming region ' num2str(iR) ' of ' num2str(iNumROI) '...'])
    
    % Find windows belonging to current region
    vWinIdx = find(vROI == vUniqueROI(iR));
    
    % Initialize
    mTi = []; mTc = [];
    
    % Loop through each window in current region
    for iW = 1:numel(vWinIdx)
        
        % Index
        iIdx = vWinIdx(iW);
        
        % If current window is not the chosen window
        if iIdx ~= iWinIdx
            
            % Estimate transformation from current window to coordinate
            % system of the chosen window
            if isempty(mTi)
                mTi = windowTransform(cL,cR,iWinIdx,iIdx);
            end
            
            % Transform triangulated points to coordinate system of the
            % chosen window
            mPtsT = mTi * cL{iIdx}.TriangulatedPoints;
            
        else
            
            % Get triangulated points
            mPtsT = cL{iIdx}.TriangulatedPoints;
            
        end
        
        % Transform triangulated points from coordinate system of the
        % chosen window to the UTM coordinate system (initial
        % transformation)
        mPtsT = sT.trans * mPtsT;
        
        % Apply transformations from alignment and optimization routines of
        % chosen window
        mPtsT = [mTa * mPtsT; ones(1,size(mPtsT,2))];
        mPtsT = transformUsingSolverVar(mPtsT,sOutput);
        
        % Apply correction for curvature of earth
        if isempty(mTc)
            [mPtsT,mTc] = curvatureCorrection(mPtsT,sT.zone,sT.hemi, ...
                'forward');
        else
            mPtsT = curvatureCorrection(mPtsT,sT.zone,sT.hemi, ...
                'forward',mTc);
        end
        
        % Save georeferencing data
        sGeo = cL{iWinIdx}.GeorefInfo;       
        sGeo.Initial.Source = cL{iWinIdx}.Properties.Source;
        if iIdx ~= iWinIdx
            sGeo.Initial.WindowTransform = mTi;
        else
            sGeo.Initial.WindowTransform = 'none';
        end
        sGeo.Initial.Triangulated2WorldTransform.curvature = mTc;
        cL{iIdx}.GeorefInfo = sGeo;
        
        % Save triangulated points after initial georeferencing
        cL{iIdx}.TriangulatedPointsGeoref = mPtsT;
        
    end
    
    catch objExc

    warning(objExc.message)
    warning(['Transformation of region ' num2str(iR) ...
        ' failed, skipping...'])
    
    end        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mTi = windowTransform(cL,cR,iWinIdx,iIdx)
    
% Select some stereo points from current window
mPts1 = cL{iIdx}.ImagePoints;
vIdx = randperm(length(mPts1),min([1000 length(mPts1)]));
mPts1 = mPts1(:,vIdx);
mPts2 = cR{iIdx}.ImagePoints;
mPts2 = mPts2(:,vIdx);

% Triangulate using matrices of the chosen window
mPtsT1 = triangulate(mPts1,mPts2, ...
cL{iWinIdx}.PoseMatrix,cR{iWinIdx}.PoseMatrix, ...
cL{iWinIdx}.IntrinsicMatrix, ...
cR{iWinIdx}.IntrinsicMatrix,false);
clear mPts1 mPts2

% Use Horn's method to compute transformation
mPtsT2 = cL{iIdx}.TriangulatedPoints;
mPtsT2 = mPtsT2(:,vIdx);
lIn = all(isfinite(mPtsT1),1) & all(isfinite(mPtsT2),1);
sTi = absor(mPtsT2(1:3,lIn),mPtsT1(1:3,lIn),'doScale',1);
mTi = sTi.M;
clear mPtsT1 mPtsT2

end
end