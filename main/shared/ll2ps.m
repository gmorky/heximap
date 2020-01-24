function mPts = ll2ps(mPts,vParams)
%Convert geographic lat lon to projected polar stereographic coordinates
%using WGS84 ellipsoid. Equations taken from:

%Snyder, J. P., 1987; Map Projections - A Working Manual. U.S.
%Geological Survey Professional Paper 1395.

%Army, Department of, 1973; Universal Transverse Mercator Grid, U. S. Army
%Technical Manual TM 5-241-8, 64 p. Superseded by DMATM 8358.2 The
%Universal Grids: Universal Transverse Mercator (UTM) and Universal Polar
%Stereographic (UPS).

%Get projection parameters
dLon0 = vParams(1);
dLat0 = vParams(2);
dFE = vParams(3);
dFN = vParams(4);
dK0 = vParams(5);

%Choose hemisphere
lN = mPts(:,2) > 0;
if all(lN)
    strH = 'N';
elseif ~all(lN)
    strH = 'S';
else
    error(['Input coordinates cannot contain both positive and ' ...
        'negative latitudes.'])
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
    
    %Get latitude and longitude coordinates
    vLon = mPts(vBlocks(i):vBlocks(i+1)-1,1);
    vLat = mPts(vBlocks(i):vBlocks(i+1)-1,2);
    
    %Latitude in polar stereographic is always treated as positive,
    %regardless of hemisphere
    vLat = abs(vLat);
    dLat0 = abs(dLat0);
    
    %Convert input to radians
    vLon = vLon*pi/180;
    vLat = vLat*pi/180;
    dLat0 = dLat0*pi/180;
    dLon0 = dLon0*pi/180;
    
    %Define constants
    dA = 6378137;
    dF = 1/298.257223563;
    dB = dA*(1-dF);
    dE = sqrt((dA^2-dB^2)/dA^2);
    
    %Compute parameter t
    vT = tan(pi/4-vLat/2)./((1-dE*sin(vLat))./(1+dE*sin(vLat))).^(dE/2);
    clear vLat
    
    %Compute parameter rho
    if abs(dLat0) == pi/2
        vRho = 2*dA*dK0*vT/sqrt((1+dE)^(1+dE)*(1-dE)^(1-dE));
    else
        dT0 = tan(pi/4-dLat0/2)/((1-dE*sin(dLat0))/(1+dE*sin(dLat0))) ...
            ^(dE/2);
        dM0 = cos(dLat0)/sqrt(1-dE^2*sin(dLat0)^2);
        vRho = dA*dM0*vT/dT0;
    end
    clear vT

    %Compute northing based on hemisphere
    if strcmpi(strH,'N')
        vY = dFN - vRho.*cos(vLon-dLon0);
    else
        vY = dFN + vRho.*cos(vLon-dLon0);
    end

    %Compute easting
    vX = dFE + vRho.*sin(vLon-dLon0);
    clear vLon vRho

    %Assign Output
    mPts(vBlocks(i):vBlocks(i+1)-1,1:2) = [vX vY];
    clear vX vY
    
end
