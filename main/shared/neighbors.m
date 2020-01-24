function mNeighbors = neighbors(mI,lInterp)

if lInterp
    
    % Add buffer around edges using linear interpolation
    vL = mI(:,1)-diff(mI(:,1:2),1,2);
    vR = mI(:,end)+diff(mI(:,end-1:end),1,2);
    mI = [vL mI vR];
    vT = mI(1,:)-diff(mI(1:2,:),1,1);
    vB = mI(end,:)+diff(mI(end-1:end,:),1,1);
    mI = [vT; mI; vB];

else
    
    % Add NaN buffer around edges
    mI = padarray(mI,[1 1],NaN);
    
end

% Define index matrix
vSz = size(mI);
[mC,mR] = meshgrid(2:vSz(2)-1,2:vSz(1)-1);
mIdxN = [-1 0;-1 1;0 1;1 1;1 0;1 -1;0 -1; -1 -1];

% Initialize
mNeighbors = zeros((vSz(1)-2)*(vSz(2)-2),8);

% Loop for each neighbor
for iN = 1:8
    
    % Define index
    mRn = mR + mIdxN(iN,1);
    mCn = mC + mIdxN(iN,2);
    vIdx = sub2ind(vSz,mRn(:),mCn(:));
    
    % Save neighbor pixel values
    mNeighbors(:,iN) = mI(vIdx);
    
end
