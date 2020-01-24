function cWindow = extChooseWindows(cM,hW)
% User selects regions for processing

% Update waitbar
try
waitbar(4/8,hW,'user selecting region(s) of interest...')
catch
end

% Get data from mat files
objML = cM{1};
objMR = cM{2};
mH1 = objML.LeftHomography;
mH2 = objML.RightHomography;
vSzL = fliplr(size(objML,'Image'));
vSzR = fliplr(size(objMR,'Image'));
vSzLs = fliplr(size(objML,'Image10'));
dScale = 10;iR = 1;iW = 1;cWindow = {};
strRepeat = 'yes';

% Determine image overlap
mOver = mH1 \ mH2 * [1 1 1; 1 vSzR(2) 1]';
mOver = mOver ./ repmat(mOver(3,:),3,1);
mOver = mOver(1:2,:)' / dScale;
mOver(:,1) = floor(min(mOver(:,1)));
mOver(3) = 1; mOver(4) = vSzLs(2);

% Show image
f = figure;
warning('off','images:initSize:adjustingMag');
imshow(adapthisteq(objML.Image10)),colormap(bone(256)),hold on
plot(mOver(:,1),mOver(:,2),'--','Color',[1 0.5 0])
plot([vSzLs(1) vSzLs(1)],[1 vSzLs(2)],'--','Color',[1 0.5 0])
text(vSzLs(2)/2,1,['Drag box to select region of interest. '  ...
                   'Double click inside box when finished'], ...
                   'Color',[1 0.5 0], ...
                   'Background','k', ...
                   'VerticalAlignment','top');
warning('on','images:initSize:adjustingMag');

% Initialize
iWinSzPix = 4400;
iBuffPix = 300;

while strcmpi(strRepeat,'yes') == 1

    % User moves rectangle to choose region of interest
    dWinSz = round(iWinSzPix/dScale);
    vPrevPos = [mOver(1) mOver(3) dWinSz dWinSz];
    h = imrect(gca,vPrevPos); zoom out
    setPositionConstraintFcn(h,@(x) ...
        constrainFunction(x,mOver,vSzLs,dWinSz));
    mROI = wait(h);
    
    % Plot the window and delete the graphic object
    rectangle('Position',mROI,'EdgeColor',[1 0.5 0])
    delete(h)
    
    % Define ROI (for left image)
    mROI = [mROI(1:2); mROI(3:4)] * dScale;
    mROI(2) = mROI(1) + mROI(2) - 1;
    mROI(4) = mROI(3) + mROI(4) - 1;
    mROI = round(mROI);
    
    % Define window vectors
    vRoiSz = fliplr(diff(mROI))+1;
    vX = round(linspace(mROI(1),mROI(2),round(vRoiSz(2)/iWinSzPix)+1));
    vY = round(linspace(mROI(3),mROI(4),round(vRoiSz(1)/iWinSzPix)+1));
    
    % Loop through each window
    for i = 1:length(vX)-1
        for j = 1:length(vY)-1
            
            % Define window, include edge buffer
            mWinL = [vX(i) vY(j); vX(i+1) vY(j+1)] + ...
                [-iBuffPix -iBuffPix; iBuffPix iBuffPix];
    
            % Make sure left image window boundaries are valid
            mWinL(mWinL < 1) = 1;
            mWinL(mWinL(:,1) > vSzL(1),1) = vSzL(1);
            mWinL(mWinL(:,2) > vSzL(2),2) = vSzL(2);
            
            % Transform window to right image
            mWinR = mH2 \ mH1 * [mWinL [1;1]]';
            mWinR = mWinR ./ repmat(mWinR(3,:),3,1);
            mWinR = mWinR(1:2,:)';
            mWinR = round(mWinR);
            
            % Make sure right image window boundaries are valid
            mWinR(mWinR < 1) = 1;
            mWinR(mWinR(:,1) > vSzR(1),1) = vSzR(1);
            mWinR(mWinR(:,2) > vSzR(2),2) = vSzR(2);
            
            % Make both windows the same size
            vWinSz = min([diff(mWinL);diff(mWinR)])-1;
            vCenL = mean(mWinL);
            vCenR = mean(mWinR);
            mWinL = round([vCenL-vWinSz/2;vCenL+vWinSz/2]);
            mWinR = round([vCenR-vWinSz/2;vCenR+vWinSz/2]);           
            
            % Save data in structure and store in cell array
            sWindow = struct('left',mWinL,'right',mWinR, ...
                'region',iR);
            cWindow{iW} = sWindow;
            
            % Increment window count
            iW = iW + 1;
            
        end
    end
    
    % User decides if they want to place another window
    strRepeat = questdlg(['Do you want to specify another region of ' ...
                          'interest?']);
    if isempty(strRepeat)
        strRepeat = 'yes';
    end
    
    % Increment ROI count
    iR = iR + 1;

end

% Close the figure window
close(f),pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vPos = constrainFunction(vPos,mOver,vSzLs,dWinSz)
    
    % Define boundaries
    vXlim = [mOver(2) vSzLs(1)];
    vYlim = get(gca,'YLim');
    
    % Constrain rectangle to specified pixel increments
    vPos(3) = max([round(vPos(3)/dWinSz)*dWinSz dWinSz]);
    vPos(4) = max([round(vPos(4)/dWinSz)*dWinSz dWinSz]);
    
    % Constrain rectangle to boundaries
    if vPos(1) < vXlim(1) || vPos(2) < vYlim(1) || ...
       vPos(1)+vPos(3) > vXlim(2) || vPos(2)+vPos(4) > vYlim(2)
        vPos = vPrevPos;
    end

    % Save position
    vPrevPos = vPos;
    
end
end
