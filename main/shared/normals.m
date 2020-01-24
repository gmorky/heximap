function varargout = normals(mX,mY,mZ,lInterp)

% Spacing
mXc = diff(mX,1,2); mYc = diff(mY,1,1);
dX = mXc(1); dY = mYc(1);

% Make sure grid is equally spaced
mXc = abs(mXc); mYc = abs(mYc);
dTol = mean([mXc(:);mYc(:)])*1E-9;
if ~all(mXc(:)-mXc(1)<dTol) || ~all(mYc(:)-mYc(1)<dTol)
    error('Grid must be equally spaced.')
end

% Pixel neighborhoods
mNeighbors = neighbors(mZ,lInterp);

% Gradients
vSz = size(mZ);
mP = reshape((mNeighbors(:,3)-mNeighbors(:,7))/(2*dX),vSz);
mQ = reshape((mNeighbors(:,5)-mNeighbors(:,1))/(2*dY),vSz);

% Normals
mNormals = repmat(1./(1 + mP(:).^2 + mQ(:).^2).^0.5,1,3) .* ...
    [-mP(:) -mQ(:) ones(prod(vSz),1)];

% Output
varargout{1} = mNormals;
varargout{2} = mP;
varargout{3} = mQ;
