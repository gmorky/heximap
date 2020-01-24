function [mCorners,cCornerFigs] = stiGetCorners(strPath,strFile,strH,hW)
% Detect corners of scanned Hexagon filmstrip

% Update waitbar
try
waitbar(1/6,hW,'detecting image corners...')
catch
end

% Make cell array to hold images of detected corners
cCornerFigs = cell(1,2);

% Define full and sub image sizes
sInfo = imfinfo(strcat([strPath strFile]));
iH = sInfo.Height;
iW = sInfo.Width;
iSz = 3000;

% Make sure correct structure is being accessed
if numel(sInfo) > 1
    sInfo = sInfo(1);
end

% Upper left corner of "a" image
if strcmpi(strH,'a')

    % Read image
    mI = imread(sInfo.Filename,'PixelRegion',{[1 iSz] [1 iSz]});
    
    % Make logical matrix defining shape of upper left corner
    vSz = [999 999];
    vC = floor((vSz+1)/2);
    lIn = false(vSz);
    lIn(vC(1):end,vC(2)) = true;
    lIn(vC(1),vC(2):end) = true;
    
    % Find the corner
    vCornU = findCorner(mI,lIn,iSz,1);
    
end

% Lower left corner of "a" image
if strcmpi(strH,'a')

    % Read image
    mI = imread(sInfo.Filename,'PixelRegion',{[iH-iSz+1 iH] [1 iSz]});
    
    % Make logical matrix defining shape of lower left corner
    vSz = [999 999];
    vC = floor((vSz+1)/2);
    lIn = false(vSz);
    lIn(1:vC(1),vC(2)) = true;
    lIn(vC(1),vC(2):end) = true;
    
    % Find the corner
    vCornL = findCorner(mI,lIn,iSz,2);
    vCornL(2) = vCornL(2) + (iH-iSz);
    vCornL(1) = iW;
    
end

% Upper right corner of "b" image
if strcmpi(strH,'b')

    % Read image
    mI = imread(sInfo.Filename,'PixelRegion',{[1 iSz] [iW-iSz+1 iW]});
    
    % Make logical matrix defining shape of upper right corner
    vSz = [999 999];
    vC = floor((vSz+1)/2);
    lIn = false(vSz);
    lIn(vC(1):end,vC(2)) = true;
    lIn(vC(1),1:vC(2)) = true;
    
    % Find the corner
    vCornU = findCorner(mI,lIn,iSz,1);
    vCornU(1) = 1;
    
end

% Lower right corner of "b" image 
if strcmpi(strH,'b')
    
    % Read image
    mI = imread(sInfo.Filename,'PixelRegion',{[iH-iSz+1 iH] ...
                                              [iW-iSz+1 iW]});
    
    % Make logical matrix defining shape of lower right corner
    vSz = [999 999];
    vC = floor((vSz+1)/2);
    lIn = false(vSz);
    lIn(1:vC(1),vC(2)) = true;
    lIn(vC(1),1:vC(2)) = true;
    
    % Find the corner
    vCornL = findCorner(mI,lIn,iSz,2);
    vCornL(2) = vCornL(2) + (iH-iSz);
    vCornL(1) = vCornL(1) + (iW-iSz);
    
end

% Save output
mCorners = floor([vCornU; vCornL])+1;
mCorners = mCorners + [10 10; -10 -10];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vCorner = findCorner(mI,lIn,iSz,iCN)

try
    
    % Save copy of image
    mIc = mI;

    % Set threshold to remove white text around borders of image
    mI(mI > 250) = 0;

    % Compute row and column means
    mH = repmat(blockproc(mI,[1 iSz],@(x) mean(x.data)),1,iSz);
    mV = repmat(blockproc(double(mI),[iSz 1],@(x) mean(x.data)),iSz,1);

    % Combine into single image
    mI = zeros(size(mH));
    mI(mH <= mV) = mH(mH <= mV);
    mI(mV <= mH) = mV(mV <= mH);

    % Compute threshold with Otsu's method and fill holes
    mI = imquantize(mI,multithresh(mI,1)) > 1;
    mI = imfill(~imfill(~mI,'holes'),'holes');

    % Find corners using Harris detector
    mC = corner(mI);

    % Compute xy image gradients to highlight edges
    [mGx,mGy] = imgradientxy(mI);
    mG = double(mGx ~= 0);
    mG(mGy ~= 0) = -1;

    % Compute local standard devation around each corner pixel using
    % logical matrix, then get location with maximum standard deviation
    vStd = localSTD(mG,lIn,mC);
    vCorner = mC(vStd == max(vStd),:);
    vCorner = vCorner(1,:);

    % Save figure of corner
    f = figure;
    vP = get(f,'Position');
    vP(2) = vP(2)/2;
    set(f,'Position',vP,'name',strFile,'numbertitle','off');
    imagesc(mIc),colormap(bone(256)),axis off,hold on
    scatter(vCorner(1),vCorner(2),[],[1 0.5 0],'s'),hold off
    objF = getframe(gcf);
    cCornerFigs{iCN} = objF.cdata;

    % Let user decide if corner is correct
    strC = questdlg('Is the detected corner correct?','','No');
    if strcmpi(strC,'Cancel')
        strC = 'No';
    end
    close(f)

    % Let user manually position corner if needed
    if strcmpi(strC,'No')
        vCorner = getCornerManual(mIc,iCN);
    end

catch
    
    % Let user manually position corner if an error occurs
    vCorner = getCornerManual(mIc,iCN);
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vStd = localSTD(mI,lIn,mPts)
% Compute local standard deviation at specific pixel locations in image

% Get size of logical matrix and number of points
vSzL = size(lIn);
iNumPts = size(mPts,1);

% Pad image and fix location coordinates for padded image
mI = padarray(mI,vSzL,'symmetric');
mPts = mPts + repmat(vSzL,iNumPts,1);

% Compute index variables
vCen = floor((vSzL+1)/2);
dH = vCen(1)-1; dW = vCen(2)-1;

% Loop for each location
vStd = zeros(1,iNumPts);
for i = 1:iNumPts
    
    % Get sub-image
    vIdx = fliplr(mPts(i,:));
    mIs = mI(vIdx(1)-dH:vIdx(1)+dH,vIdx(2)-dW:vIdx(2)+dW);
    
    % Compute standard deviation
    vStd(i) = std(mIs(lIn));
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vCorner = getCornerManual(mIc,iCN)
    
    % Show figure
    f = figure;
    vP = get(f,'Position');
    vP(2) = vP(2)/2;
    set(f,'Position',vP,'name',strFile,'numbertitle','off');
    imagesc(mIc),colormap(bone(256)),axis off,hold on
    text(size(mIc,2)/2,1, ...
        'Click on image corner: ', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','Top', ...
        'Color',[1 0.5 0], ...
        'Background','k')
    
    % User inputs corner point
    [iX,iY,~] = ginput(1);
    vCorner = [iX iY];
    scatter(vCorner(1),vCorner(2),[],[1 0.5 0],'s'),hold off
    
    % Save figure of corner
    objF = getframe(gcf);
    cCornerFigs{iCN} = objF.cdata;
    close(f)
    
end
end