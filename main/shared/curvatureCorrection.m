function [mPtsT,varargout] = curvatureCorrection(mPtsT,iZ,strH,strDir, ...
    varargin)
% When converting from the coordinate system of the triangulated points
% directly to UTM (using a rigid 3D transformation), vertical errors are
% introduced due to the curvature of the earth. This function corrects the
% errors using an intermediate transformation to the ECEF coordinate
% system.

if strcmpi(strDir,'forward')

    % Choose random subset of triangulated points (they are in UTM)
    mPtsTs = geoSamplePoints(mPtsT,100);

    % Transform points from UTM to WGS84 to ECEF
    mPtsWs = utm2ll(mPtsTs',iZ,strH)';
    [mPtsWs(1,:),mPtsWs(2,:),mPtsWs(3,:)] = geodetic2ecef( ...
        referenceEllipsoid('wgs84'),mPtsWs(2,:),mPtsWs(1,:),mPtsWs(3,:));

    % Rigid transformation from UTM to ECEF
    try
        mT = varargin{1};
    catch
        lInC = all(isfinite(mPtsTs),1) & all(isfinite(mPtsWs),1);
        sTc = absor(mPtsTs(1:3,lInC),mPtsWs(1:3,lInC),'doScale',1);
        mT = sTc.M;
    end

    % Apply the rigid transformation
    mPtsT = mT * mPtsT;

    % Transform points from ECEF back to WGS84 back to UTM (non-rigid)
    [mPtsT(2,:),mPtsT(1,:),mPtsT(3,:)] = ecef2geodetic( ...
        referenceEllipsoid('wgs84'),mPtsT(1,:),mPtsT(2,:),mPtsT(3,:));
    mPtsT = ll2utm(mPtsT',iZ,strH)';

elseif strcmpi(strDir,'inverse')
    
    % Transform points from UTM to WGS84 to ECEF
    mPtsT = utm2ll(mPtsT',iZ,strH)';
    [mPtsT(1,:),mPtsT(2,:),mPtsT(3,:)] = geodetic2ecef( ...
    referenceEllipsoid('wgs84'),mPtsT(2,:),mPtsT(1,:),mPtsT(3,:));

    % Undo rigid transformation (from ECEF back to UTM)
    try
        mT = varargin{1};
    catch
        error('Curvature transform matrix required for inverse transform.')
    end
    mPtsT = mT \ mPtsT;
    
else
    error('Invalid transformation direction parameter.')
end

% Specify rigid transform matrix as optional output
varargout{1} = mT;
