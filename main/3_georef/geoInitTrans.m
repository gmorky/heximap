function [] = geoInitTrans(objL,strRef,strShpPath,lVis,hW)
% Compute transformations for initial georeferencing (using chosen window)

% Update waitbar
try
waitbar(5/8,hW,'computing transformations for initial georeferencing...')
catch
end

% Update command window
disp(['computing transformations for initial georeferencing using' ...
     ' window ' num2str(objL.WindowID) '...'])

% Get ground control points
sGeo = objL.GeorefInfo;
mPtsT = sGeo.Initial.GroundControlPoints.Triangulated;
mPtsW = sGeo.Initial.GroundControlPoints.World;

% Convert world points to UTM
[mPtsW,iZ,strH] = ll2utm(mPtsW,[],[]);

% Use Horn's absolute orientation method to compute initial transformation
mPtsW = mPtsW'; mPtsT = mPtsT';
lIn = all(isfinite(mPtsT),1) & all(isfinite(mPtsW),1);
sT = absor(mPtsT(1:3,lIn),mPtsW(1:3,lIn),'doScale',1);
mT = sT.M;

% Choose random subset of triangulated points (to conserve memory)
mPtsT = geoSamplePoints(mT*objL.TriangulatedPoints,1E6);

% Improve accuracy using nonlinear optimization
sOpt.rotation = [45 45 45];
sOpt.translation = [10E3 10E3 10E3];
sOpt.scale = [NaN NaN NaN];
sOpt.globalScale = 2;
sOpt.maxIterations = 200;
sOpt.visualize = lVis;
sOpt.polySurf = false;
[mTa,sOutput] = geoOptiTrans(mPtsT,strRef,strShpPath,iZ,strH,sOpt);

% Make sure vertical RMSE error is reasonable
if sOutput.verticalRMSE > 50
    error('The vertical RMSE is high. Try again with a different window.')
end

% Save in mat file
sT = struct();
sT.trans = mT;
sT.zone  = iZ;
sT.hemi = strH;
sGeo = objL.GeorefInfo;
sGeo.Initial.Triangulated2WorldTransform = sT;
sGeo.Initial.AlignmentOutput = mTa;
sGeo.Initial.OptimizationOutput = sOutput;
objL.GeorefInfo = sGeo;
