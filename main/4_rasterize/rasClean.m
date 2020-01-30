function mDem = rasClean(mDem,dResM,lClean,iGap,iSpec)

if lClean

    % Interpolate small data gaps
    iGap = round(iGap/dResM^2);
    lNaN = isnan(mDem);
    mDem = inpaint_nans(mDem);
    lKeep = bwareaopen(lNaN,iGap);
    mDem(lKeep) = NaN;

    % Dilate remaining data gaps
    lKeep = imdilate(lKeep,strel('disk',1));
    mDem(lKeep) = NaN;

    % Remove isolated speckles
    iSpec = round(iSpec/dResM^2);
    mDem(~bwareaopen(~isnan(mDem),iSpec)) = NaN;

end
