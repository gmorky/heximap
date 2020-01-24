function [mPts,varargout] = ll2utm(mPts,iZ,strH)
%Convert geographic lat lon to projected UTM coordinates using WGS84
%ellipsoid. Equations taken from:

%Snyder, J. P., 1987; Map Projections - A Working Manual. U.S.
%Geological Survey Professional Paper 1395.

%Army, Department of, 1973; Universal Transverse Mercator Grid, U. S. Army
%Technical Manual TM 5-241-8, 64 p. Superseded by DMATM 8358.2 The
%Universal Grids: Universal Transverse Mercator (UTM) and Universal Polar
%Stereographic (UPS).

%Also see http://www.uwgb.edu/dutchs/usefuldata/utmformulas.htm

% Find non-NaN points
lIn = ~any(isnan(mPts),2);

%Choose zone if necessary
if isempty(iZ)
    dLon0 = floor(mean(mPts(lIn,1))/6)*6+3;
    iZ = floor(dLon0/6)+31;
end

%Choose hemisphere if necessary
if isempty(strH)
    if mean(mPts(lIn,2)) > 0
        strH = 'N';
    else
        strH = 'S';
    end
end

%Assign false northing
if strH == 'N';
    dFN = 0;
elseif strH == 'S';
    dFN = 1E7;
else
    error('Hemisphere must be N or S.')
end

%Compute intervals for block processing to conserve memory
iSz = size(mPts,1);
iNumWinX = ceil(iSz/5E6);
vBlocks = round(linspace(1,iSz,iNumWinX));
vBlocks(end) = vBlocks(end)+1;

%Fix if number of points is smaller than one block
if length(vBlocks) == 1
    vBlocks = [1 vBlocks];
end

%Loop for each block
for i = 1:length(vBlocks)-1
    
    %Get latitude and longitude coordinates in radians
    vLon = mPts(vBlocks(i):vBlocks(i+1)-1,1) * pi/180;
    vLat = mPts(vBlocks(i):vBlocks(i+1)-1,2) * pi/180;

    %Define constants
    dA = 6378137;
    dB = 6356752.314245;
    dK0 = 0.9996;
    dE = sqrt(1-dB^2/dA^2);
    dE2 = (dE*dA/dB)^2;
    dLon0 = (6.*iZ-183) * pi/180;
    
    %Compute meridional arc
    vM = dA*((1 - dE^2/4 - 3*dE^4/64 - 5*dE^6/256) * vLat ...
        - (3*dE^2/8 + 3*dE^4/32 + 45*dE^6/1024) * sin(2*vLat) ... 
        + (15*dE^4/256 + 45*dE^6/1024) * sin(4*vLat) ...
        - (35*dE^6/3072) * sin(6*vLat));
    
    %Define some more variables
    vN = dA./sqrt(1-dE.^2*sin(vLat).^2); 
    vT = tan(vLat).^2;                
    vC = (dE.^2)/(1-dE.^2).*cos(vLat).^2;
    vA = (vLon-dLon0).*cos(vLat);
    clear vLon

    %Convert longitude to easting
    vX = dK0.*vN.*(vA+(1-vT+vC).*vA.^3./6+(5-18.*vT+vT.^2+72.*vC-58.*dE2).* ...
        vA.^5./120);
    
    %Convert latitude to northing
    vY = dK0.*vM+dK0.*vN.*tan(vLat).*(vA.^2./2+(5-vT+9.*vC+4.*vC.^2).* ...
        vA.^4./24+(61-58*vT+vT.^2+600*vC-330.*dE2).*vA.^6./720);
    clear vM vN vT vC vA vLat
    
    %Assign Output
    mPts(vBlocks(i):vBlocks(i+1)-1,1:2) = [vX+5E5 vY+dFN];
    clear vX vY

end

%Output zone and hemisphere
varargout{1} = iZ;
varargout{2} = strH;
