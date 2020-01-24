function [] = extTriangulate(objL,objR,hW,cReg)
%Triangulate points in 3D space, use block processing to conserve memory

%Update waitbar
try
set(get(findobj(hW,'type','axes'),'title'), 'string', ...
    ['region ' cReg{1} ' of ' cReg{2} ': triangulating points...'])
pause(0.1)
catch
end

%Define blocks
vSz = size(objL,'ImagePoints');
iWinSz = 1E6;
iNumWin = ceil(vSz(2)/iWinSz);
if iNumWin < 2
    iNumWin = 2;
end
vB = round(linspace(1,vSz(2),iNumWin));

%Get camera data
mP1 = objL.PoseMatrix;
mP2 = objR.PoseMatrix;
mK1 = objL.IntrinsicMatrix;
mK2 = objR.IntrinsicMatrix;

%Loop for each block
objL.TriangulatedPoints = [];
for iB = 1:iNumWin-1
    
    %Triangulate the points
    objL.TriangulatedPoints(1:4,vB(iB):vB(iB+1)) = triangulate( ...
        objL.ImagePoints(1:3,vB(iB):vB(iB+1)), ...
        objR.ImagePoints(1:3,vB(iB):vB(iB+1)), ...
        mP1,mP2,mK1,mK2,false);
    
    % Display progress
    disp(['  ' num2str(round(iB/(iNumWin-1)*100)) '% complete...'])

end
