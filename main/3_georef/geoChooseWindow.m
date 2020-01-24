function iWinIdx = geoChooseWindow(cL,cR,hW)
% User selects window for initial georeferencing

% Update waitbar
try
waitbar(3/8,hW,'user selecting window for initial georeferencing...')
catch
end

% Make sure the windows were all extracted from the same source images
cInfoL = cellfun(@(x) x.SourceImageInfo,cL,'Uni',0);
cInfoR = cellfun(@(x) x.SourceImageInfo,cR,'Uni',0);
if ~isequal(cInfoL,repmat(cInfoL(1),size(cInfoL))) || ...
   ~isequal(cInfoR,repmat(cInfoR(1),size(cInfoR)))
    error(['Some files in the selected folder were not extracted ' ...
        'from the same source image pair.'])
end

% Show image
hF1 = figure;
warning('off','images:initSize:adjustingMag');
imshow(adapthisteq(cL{1}.SourceImage)),colormap(bone(256)),hold on
vSz = fliplr(size(cL{1},'SourceImage'));
text(vSz(2)/2,1,['The orange boxes indicate regions where '  ...
                   'DEMs have been extracted.'], ...
                   'Color',[1 0.5 0], ...
                   'Background','k', ...
                   'VerticalAlignment','top');
warning('on','images:initSize:adjustingMag');

% Get labels
vLabels = cell2mat(cellfun(@(x) x.WindowID,cL,'Uni',0));

% Plot the windows and their labels
for i = 1:numel(cL)
    
    % Rectangle
    dScale = 10;
    mWin = cL{i}.Window/dScale;
    vRec = [min(mWin(:,1)) min(mWin(:,2)) ...
        abs(diff(mWin(:,1))) abs(diff(mWin(:,2)))];
    rectangle('Position',vRec,'EdgeColor',[1 0.5 0])
    
    % Label
    text(mean(mWin(:,1)),mean(mWin(:,2)),num2str(vLabels(i)), ...
        'Color',[1 0.5 0], ...
        'Background','k', ...
        'VerticalAlignment','middle')
    
end

% Initialize
lValidInput = false(1,2);

% Make sure user input is valid
while any(~lValidInput)
    
    % Input dialog box
    strP = 'Which processing window to use for initial georeferencing?';
    strTitle = 'Input';
    iLabel = str2double(inputdlg(strP,strTitle,1));
    iWinIdx = find(iLabel == vLabels);
    
    % If user hit cancel button
    if isempty(iLabel)
        error('Process aborted by user.')
    end
    
    % Make sure window exists, and no previous errors occured in the chosen
    % window (i.e. it has the data needed.)
    if sum(iLabel == vLabels) == 1
        lValidInput(1) = true;
        if ~isempty(whos(cL{iWinIdx},'TriangulatedPoints'))
            lValidInput(2) = true;
        end        
    end

    % Exit while-loop if input is valid
    if all(lValidInput)
        continue
    end
    
    % Warn
    uiwait(warndlg(['Input is invalid, or errors occurred during ' ...
        'previous extraction step for the chosen window. ' ...
        'Please try again.'],'Warning'));
    lValidInput = false(1,2);
    
end

% Read the image and disparity map, resample to save memory
mI = imresize(cL{iWinIdx}.RectImage,0.5,'bilinear');
mD = imresize(cL{iWinIdx}.Disparity,size(mI),'nearest');

% Find disparity map holes
lS = repmat(isnan(mD),[1 1 3]);

% Make image combined with transparent disparity map
warning('off','images:initSize:adjustingMag')
figure('visible','off')
mI = im2uint8(ind2rgb(gray2ind(mat2gray(mI),256),bone(256)));
mD = im2uint8(ind2rgb(gray2ind(mat2gray(mD),256),jet(256)));
mIf = imfuse(mI,mD,'blend');
mIf(lS) = mI(lS);
clear lS mD

% Show figure
hF2 = figure;
hI = imshow(mIf);
colormap(bone(256))
warning('on','images:initSize:adjustingMag')
title(['Region ' num2str(iLabel) ' disparity map (click to toggle)'])

% User input to toggle disparity map
lHidden = false;
set(hI,'ButtonDownFcn',{@toggleDisparity,mI,mIf})

% Show message
warndlg(['Make sure this region you selected has some vertical relief, ' ...
    'a good disparity map, and doesn''t contain a significant amount of ' ...
    'unstable terrain.'], ...
    'Warning');

% Create "use this region" push button
iW = 0.2; iH = 0.05; iL = (0.5-iW/2)-iW; iB = 0.05;
uicontrol(hF2,'Style','pushbutton','String','Use this region',...
    'Units','normalized', ...
    'Position',[iL iB iW iH],...
    'Callback',{@thisRegionCallback,hF2});

% Create "choose new region" push button
iW = 0.2; iH = 0.05; iL = (0.5-iW/2)+iW; iB = 0.05;
uicontrol(hF2,'Style','pushbutton','String','Choose new region',...
    'Units','normalized', ...
    'Position',[iL iB iW iH],...
    'Callback',{@newRegionCallback,hF2});

% Wait for user to click on a push button (which closes the figure)
uiwait(hF2)

% Close the first figure window
close(hF1),pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = toggleDisparity(src,~,mI,mIf)
    
    % Set image data
    if lHidden
        set(src,'CData',mIf)
    else
        set(src,'CData',mI)
    end
    
    % Toogle flag
    lHidden = ~lHidden;
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = thisRegionCallback(~,~,hF)
    close(hF)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = newRegionCallback(~,~,hF)
    close(hF)
    iWinIdx = '';
end
end
