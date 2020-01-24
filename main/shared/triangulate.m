function varargout = triangulate(mPts1,mPts2,mP1,mP2,mK1,mK2,lE)
% Triangulate points in 3D space using direct linear method (Hartley and
% Zisserman)

% Save copy of original points
if lE
    mPts1O = mPts1;
    mPts2O = mPts2;
end

% Compute 3D image points
mPts1 = mK1 \ mPts1;
mPts2 = mK2 \ mPts2;

% Initialize
mPts3D = zeros(4,size(mPts1,2));

% Loop for each image point
for i = 1:size(mPts1,2)

    % Define matrix A
    mA = [mPts1(1,i) * mP1(3,:) - mP1(1,:);
          mPts1(2,i) * mP1(3,:) - mP1(2,:);
          mPts2(1,i) * mP2(3,:) - mP2(1,:);
          mPts2(2,i) * mP2(3,:) - mP2(2,:)];

    % Normalize rows
    mAn = [mA(1,:) / norm(mA(1,:));
           mA(2,:) / norm(mA(2,:));
           mA(3,:) / norm(mA(3,:));
           mA(4,:) / norm(mA(4,:))];

    % Compute 3D point as the unit singular vector corresponding to the
    % smallest singular value of matrix A
    [~,~,mV] = svd(mAn);
    mPts3D(:,i) = mV(:,end);
    
    % Normalize the point
    mPts3D(:,i) = mPts3D(:,i) / mPts3D(4,i);

end

if lE
    
    % Project the 3D points back onto image plane 1
    mPts1R = mK1 * mP1(1:3,:) * mPts3D;
    mPts1R = mPts1R ./ repmat(mPts1R(3,:),3,1);

    % Project the 3D points back onto image plane 2
    mPts2R = mK2 * mP2(1:3,:) * mPts3D;
    mPts2R = mPts2R ./ repmat(mPts2R(3,:),3,1);

    % Compute distances
    vE1 = sqrt((mPts1O - mPts1R).^2);
    vE2 = sqrt((mPts2O - mPts2R).^2);
    vError = [vE1(1,:) vE1(2,:) vE2(1,:) vE2(2,:)];
    vError(isnan(vError)) = 0;
    
else
    
    vError = [];

end

% Output
varargout{1} = mPts3D;
varargout{2} = vError;
