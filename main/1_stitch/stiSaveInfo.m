function [] = stiSaveInfo(objM,cCornerFigsA,cCornerFigsB)
% Save image info in mat file
   
% Save camera scan resolution, focal length, and estimated principal
% point
iScanRes = 7;
objM.ScanResMicrometers = iScanRes;
objM.FocalLengthMicrometers = 304.8E3;
objM.FocalLengthPixels = 304.8E3 / iScanRes;
objM.PrincipalPointPixels = fliplr(size(objM,'Image')/2);

% Save corner figures (to see if corners were detected correctly)
objM.CornerFigs = [cCornerFigsA{1} cCornerFigsB{1}; ...
                   cCornerFigsA{2} cCornerFigsB{2}];
