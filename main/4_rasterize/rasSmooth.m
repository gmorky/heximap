function mDem = rasSmooth(vX,vY,mDem,lMed,lDen,iMed,iDenT,iDenN)

if lMed
    
    % Interpolate data gaps
    lNaN = isnan(mDem);
    mDem = inpaint_nans(mDem);

    % Median filter
    mDem = medfilt2(mDem,[iMed iMed],'symmetric');
    
    % Reset data gaps
    mDem(lNaN) = NaN;

end

if lDen

    % Mesh denoising algorithm
    sInput.lon = vX;
    sInput.lat = vY;
    sInput.null = 'dem';
    sInput.params = [iDenT iDenN iDenN*2];
    sInput.blockSize = 1000;
    mDem = meshDenoise(mDem,sInput);

end
