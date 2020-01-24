function mPts = ps2ll(mPts,vParams)
%Convert projected polar stereographic to geographic lat lon coordinates
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
if dLat0 >= 0
    strH = 'N';
else
    strH = 'S';
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

    %Get easting and northing coordinates
    vX = mPts(vBlocks(i):vBlocks(i+1)-1,1);
    vY = mPts(vBlocks(i):vBlocks(i+1)-1,2);
 
    %Latitude in polar stereographic is always treated as positive,
    %regardless of hemisphere
    dLat0 = abs(dLat0);
    
    %Convert input to radians
    dLat0 = dLat0*pi/180;
    dLon0 = dLon0*pi/180;
    
    %Define constants
    dA = 6378137;
    dF = 1/298.257223563;
    dB = dA*(1-dF);
    dE = sqrt((dA^2-dB^2)/dA^2);
    
    %Subtract false easting and northing
    vX = vX - dFE;
    vY = vY - dFN;

    %Compute longitude
    if strcmpi(strH,'N')
        vLon = dLon0 + atan2(vX,-vY);
    else
        vLon = dLon0 + atan2(vX,vY);
    end
    
    %Compute parameters
    dT0 = tan(pi/4-dLat0/2)/((1-dE*sin(dLat0))/(1+dE*sin(dLat0))) ...
        ^(dE/2);
    dM0 = cos(dLat0)/sqrt(1-dE^2*(sin(dLat0))^2);
    vRho = sqrt(vX.^2+vY.^2);
    clear vX vY
    
    %Compute parameter t
    if abs(dLat0) == pi/2
        vT = vRho*sqrt((1+dE)^(1+dE)*(1-dE)^(1-dE))/(2*dA*dK0);
    else
        vT = vRho*dT0/(dA*dM0);
    end
    clear vRho

    %Compute isometric latitude
    vChi = pi/2 - 2*atan(vT);
    clear vT
    
    %Constants for converting isometric to geographic latitude
    dA = dE^2/2 + 5*dE^4/24 + dE^6/12 + 13*dE^8/360;
    dB = 7*dE^4/48 + 29*dE^6/240 + 811*dE^8/11520;
    dC = 7*dE^6/120 + 81*dE^8/1120;
    dD = 4279*dE^8/161280;
    
    %Compute geographic latitude using series
    vLat = vChi+dA*sin(2*vChi)+dB*sin(4*vChi)+dC*sin(6*vChi)+dB*sin(8*vChi);
    clear vChi
     
    if strcmpi(strH,'N')
        
        %Fix longitude phasing for northern hemisphere
        %vLon = pi-mod(vLon,2*pi);
        vLon=mod(vLon+pi,2*pi)-pi;
        
    else
        
        %Make latitude negative for southern hemisphere
        vLat = -vLat;
        
    end
    
    %Convert to degrees and assign output
    mPts(vBlocks(i):vBlocks(i+1)-1,1:2) = [vLon vLat]*180/pi;
    clear vLon vLat
    
end
