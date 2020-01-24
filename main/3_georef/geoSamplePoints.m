function [mPts,varargout] = geoSamplePoints(mPts,iNumPts)
% Choose random sample of points

% Initialize
lT = false;
vSz = size(mPts);

% Check input
if any(vSz < 2)
    error(['Must specify more than 1 point, ' ...
        'and points must have at least 2 dimensions.'])
end
if vSz(1) == vSz(2)
    warning('Ambiguous input. Assuming each column is a point.')
end

% Transpose if necessary so each column is a point
if vSz(1) > vSz(2)
    mPts = mPts';
    lT = true;
end

% Get index
vIdx = randperm(length(mPts),min([iNumPts length(mPts)]));

% Choose points
mPts = mPts(:,vIdx);

% Invert transpose if necessary
if lT
    mPts = mPts';
end

% Output
varargout{1} = vIdx;
