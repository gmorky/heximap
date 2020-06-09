function [] = stiResize(objM,hW)
% Save lower resolution copies of Hexagon image

% Update waitbar
try
waitbar(5/6,hW,'saving copies of image at lower resolutions...')
catch
end

% Scales for downsampling. Note that each successive scale is relative to
% the preceding one.
vScale = [1/2 1/5];

% Field names for saving downscaled images
cFields = {'Image','Image2','Image10'};

% Loop through each scale
for i = 1:length(vScale)
    
    % Get size of higher resolution image
    [iH,iW] = size(objM,cFields{i});

    % Spatial referencing info for higher resolution image
    sR.DeltaLon = 1;
    sR.DeltaLat = -1;
    sR.Lonlim = [1 iW] + [-0.5 0.5];
    sR.Latlim = [1 iH] + [-0.5 0.5];
    sInfo.SpatialRef = sR;
    sInfo.Width = iW;
    sInfo.Height = iH;
    sInfo.matSource = objM.Properties.Source;
    sInfo.matField = cFields{i};

    % Spatial referencing vectors for downscaled image
    dS = vScale(i);  
    vSz = round(dS*[iH iW]);
    vX = 1:vSz(2);
    vY = fliplr(1:vSz(1));

    % Resampling parameters
    sParams.blockSize = round(1000*dS);  
    sParams.matSource = objM.Properties.Source;
    sParams.matField = cFields{i+1};
    sParams.transform = [dS 0 0 0; 0 dS 0 0; 0 0 1 0; 0 0 0 1];   
    
    % Initialize field with uint8 data type
    objM.(cFields{i+1}) = uint8(0);
    
    % Resample
    warning('off','block_process:class_match')
    grid2grid(sInfo,vX,vY,sParams);
    warning('on','block_process:class_match')
    
end
