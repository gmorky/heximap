function [] = extReadImage(objM,objW,hW,cWin)
% Read Hexagon image windows

% Update waitbar
try
set(get(findobj(hW,'type','axes'),'title'), 'string', ...
    ['window ' cWin{1} ' of ' cWin{2} ': reading images...'])
pause(0.1)
catch
end

% Get window
mWin = objW.Window;

% Read the Hexagon image
objW.Image = objM.Image(mWin(3):mWin(4),mWin(1):mWin(2));
