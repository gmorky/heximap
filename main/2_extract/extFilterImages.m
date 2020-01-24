function [mI1,mI2] = extFilterImages(mI1,mI2,iE,sOpt)
% Match image histograms, apply locally adaptive contrast and noise
% filters

% Create masks for empty pixels
lM1 = mI1 ~= iE;
lM2 = mI2 ~= iE;

if sOpt.histmatch

    % Match image histograms
    mI2(lM2) = imhistmatch(mI2(lM2),mI1(lM1),256);

end

if sOpt.adapthisteq
    
    % Apply locally adaptive histogram equalization
    mI1 = adapthisteq(mI1,'ClipLimit',0.03,'NumTiles',[20 20]);
    mI2 = adapthisteq(mI2,'ClipLimit',0.03,'NumTiles',[20 20]);
    
end

if sOpt.wiener2

    % Apply noise removal filter
    mI1 = wiener2(mI1,[3 3]);
    mI2 = wiener2(mI2,[3 3]);

end

% Reassign empty pixels
mI1(~lM1) = iE;
mI2(~lM2) = iE;
