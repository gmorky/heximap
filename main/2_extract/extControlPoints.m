function [] = extControlPoints(objM,strHexFile,hW)
% User inputs corner points to roughly align Hexagon images with the
% reference DEM.  Output is a spatial transformation structure

% Update waitbar
try
waitbar(1/8,hW,'user specifying corner control points...')
catch
end

% Get control points array, or make empty one if it doesn't exist
try
    mPts = objM.CornerGCPs;
catch
    mPts = [];
end
 
% User inputs corner coordinates if needed
if isempty(mPts)
    inputCornerGCPs(objM,strHexFile)
else
    
    % Ask user if they want to enter new corner GCPs
    strN = questdlg(['Do you want to enter new corner ground' ...
                     ' control points for ' strHexFile '?'], ...
                     '','No');
    if strcmpi(strN,'Yes')
        inputCornerGCPs(objM,strHexFile)
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = inputCornerGCPs(objM,strHexFile)
    
% Input corner coordinates from metadata (see earthexplorer.usgs.gov)
cPrompt = {'NW Corner Latitude (dec):', 'NW Corner Longitude (dec):', ...
           'NE Corner Latitude (dec):', 'NE Corner Longitude (dec):', ...
           'SE Corner Latitude (dec):', 'SE Corner Longitude (dec):', ...
           'SW Corner Latitude (dec):', 'SW Corner Longitude (dec):'};
cDef = cell(size(cPrompt));
for i = 1:numel(cDef)
    cDef{i} = '';
end
mSz = repmat([1 65],8,1);       
cIn = inputdlg(cPrompt,['Image ID: ' strHexFile(1:end-4)],mSz,cDef,'on');
mPtsWld = str2double([cIn(2:2:end) cIn(1:2:end)]);

% Put corner coordinates in correct order. This assumes the Hexagon
% image was scannned with left side = north, right side = south.
mPtsWld = [mPtsWld(2,:);mPtsWld(4,:);mPtsWld(3,:);mPtsWld(1,:)];

% Define image size and corners. Also note that the y-coordinates are
% made negative to correspond with right-handed world coordinate system
[iH,iW] = size(objM,'Image');
mPtsImg = [1 1; iW iH; iW 1; 1 iH];
mPtsImg(:,2) = -mPtsImg(:,2); 

% Compute spatial transformation from control point pairs. Can ignore
% condition number warning.
warning('off','images:maketform:conditionNumberofAIsHigh')
sT = cp2tform(mPtsImg,mPtsWld,'nonreflective similarity');
warning('on','images:maketform:conditionNumberofAIsHigh')

% Save spatial transformation structure and corner coordinates to mat
% file
objM.SpatialTrans = affine2d(sT.tdata.T);
objM.CornerGCPs = mPtsWld;

end
end