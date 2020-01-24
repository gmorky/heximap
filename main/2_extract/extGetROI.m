function [cL,cR] = extGetROI(cFileL,cFileR,vROI,iROI)
% Get windows belonging to an ROI

% Get windows belonging to current ROI
lIn = vROI == iROI;
cL = cFileL(lIn);
cR = cFileR(lIn);

% Only keep mat files where stereo processing was successful
lIn = ...
    cell2mat(cellfun(@(x) ~isempty(whos(x,'ImagePoints')),cL,'Uni',0)) ...
  & cell2mat(cellfun(@(x) ~isempty(whos(x,'ImagePoints')),cR,'Uni',0));
cL = cL(lIn);
cR = cR(lIn); 

% Make sure region has at least one window
if isempty(cL) || isempty(cR)
    error(['Disparity maps were not successfully computed for any ' ...
        'windows belonging to the current region.'])
end
