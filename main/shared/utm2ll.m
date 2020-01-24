function mPts = utm2ll(mPts,iZ,strH)
%Convert projected UTM coordinates to geographic lat lon using WGS84
%ellipsoid. Equations taken from:

%Snyder, J. P., 1987; Map Projections - A Working Manual. U.S.
%Geological Survey Professional Paper 1395.

%Army, Department of, 1973; Universal Transverse Mercator Grid, U. S. Army
%Technical Manual TM 5-241-8, 64 p. Superseded by DMATM 8358.2 The
%Universal Grids: Universal Transverse Mercator (UTM) and Universal Polar
%Stereographic (UPS).

%Also see http://www.uwgb.edu/dutchs/usefuldata/utmformulas.htm

%Error if user does not specify zone
if isempty(iZ)
    error('Must specify UTM zone for the input coordinates.')
end

%Error if user does not specify hemisphere
if isempty(strH)
    error('Must specify UTM hemisphere for the input coordinates.')
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

    %Get easting and northing coordinates
    vX = mPts(vBlocks(i):vBlocks(i+1)-1,1);
    vY = mPts(vBlocks(i):vBlocks(i+1)-1,2);

    %False easting
    vX = vX - 5E5;

    %Define constants
    dA = 6378137;
    dB = 6356752.3142;
    dK0 = 0.9996;
    dE = sqrt(1-dB^2/dA^2);
    dE2 = (dE*dA/dB)^2;
    dLon0 = (6.*iZ-183) * (pi/180);

    %Compute meridional arc
    vM = (vY-dFN)/dK0;

    %Compute footprint latitude
    vMu = vM/(dA*(1-dE^2/4-3*dE^4/64-5*dE^6/256));
    dE1 = (1-(1-dE^2)^0.5)/(1+(1-dE^2)^0.5);
    dJ1 = 3*dE1/2-27*dE1^3/32;
    dJ2 = 21*dE1^2/16-55*dE1^4/32;
    dJ3 = 151*dE1^3/96;
    dJ4 = 1097*dE1^4/512;
    vFP = vMu + dJ1*sin(2*vMu) + dJ2*sin(4*vMu) + ...
          dJ3*sin(6*vMu) + dJ4*sin(8*vMu);
    clear vM vMu

    %Compute parameters for latitude and longitude
    vC1 = dE2 * cos(vFP).^2;
    vT1 = tan(vFP).^2;
    vR1 = dA*(1-dE^2)./(1-dE^2*sin(vFP).^2).^1.5;
    vN1 = dA./(1-dE^2*sin(vFP).^2).^0.5;
    vD = vX./(vN1*dK0);

    %Compute latitude
    vQ1 = vN1.*tan(vFP)./vR1;
    clear vRi vN1
    vQ2 = vD.^2/2;
    vQ3 = (5+3*vT1+10*vC1-4*vC1.^2-9*dE2).*vD.^4/24;
    vQ4 = (61+90*vT1+298*vC1+45*vT1.^2-3*vC1.^2-252*dE2).*vD.^6/720;
    vLat = (vFP-vQ1.*(vQ2-vQ3+vQ4)) * (180/pi);
    clear vQ1 vQ2 vQ3 vQ4

    %Compute longitude
    vQ5 = vD;
    vQ6 = (1+2*vT1+vC1).*vD.^3/6;
    vQ7 = (5-2*vC1+28*vT1-3*vC1.^2+8*dE2+24*vT1.^2).*vD.^5/120;
    vLon = (dLon0+(vQ5-vQ6+vQ7)./cos(vFP)) * (180/pi);
    clear vQ5 vQ6 vQ7 vFP vD

    %Assign output
    mPts(vBlocks(i):vBlocks(i+1)-1,1:2) = [vLon vLat];

end
