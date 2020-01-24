function [vX,vY] = makeSpatialRefVecs(sR,strType)

try
    
    % Fix spatial referencing if raster interpretation is 'postings'
    try
        if strcmpi(sR.RasterInterpretation,'postings')
            sR = georefcells(sR.Latlim,sR.Lonlim,sR.RasterSize);
            sR.ColumnsStartFrom = 'north';
        end
    catch
    end
    
    if strcmp(strType,'full')
    
        % Full vectors
        dX = sR.DeltaLon; dY = sR.DeltaLat;
        vX = sR.Lonlim(1) + dX/2 : dX : sR.Lonlim(2) - dX/2;
        vY = sR.Latlim(2) + dY/2 : dY : sR.Latlim(1) - dY/2;
        
    elseif strcmp(strType,'limits')
        
        % Only first and last elements of vectors
        dX = sR.DeltaLon; dY = sR.DeltaLat;
        vX = [sR.Lonlim(1) + dX/2 sR.Lonlim(2) - dX/2];
        vY = [sR.Latlim(2) + dY/2 sR.Latlim(1) - dY/2];
        
    else
        error('Unknown character input.')
    end
    
catch
    
    % Fix spatial referencing if raster interpretation is 'postings'
    try
        if strcmpi(sR.RasterInterpretation,'postings')
            sR = maprefcells(sR.XWorldLimits,sR.YWorldLimits,sR.RasterSize);
            sR.ColumnsStartFrom = 'north';
        end
    catch
    end
    
    if strcmp(strType,'full')
    
        % Full vectors
        dX = sR.DeltaX; dY = sR.DeltaY;
        vX = sR.XLimWorld(1) + dX/2 : dX : sR.XLimWorld(2) - dX/2;
        vY = sR.YLimWorld(2) + dY/2 : dY : sR.YLimWorld(1) - dY/2; 

    elseif strcmp(strType,'limits')
        
        % Only first and last elements of vectors
        dX = sR.DeltaX; dY = sR.DeltaY;
        vX = [sR.XLimWorld(1) + dX/2 sR.XLimWorld(2) - dX/2];
        vY = [sR.YLimWorld(2) + dY/2 sR.YLimWorld(1) - dY/2]; 

    else
        error('Unknown character input.')
    end
 
end
