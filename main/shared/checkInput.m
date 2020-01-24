function [] = checkInput(strFileIn)
% Check input data for possible errors

% Get file extension
[~,~,strExt] = fileparts(strFileIn);

switch strExt   
    case '.tif'

        % Make sure the raster uses standard WGS84
        sInfo = geotiffinfo(strFileIn);
        if ~strcmpi(sInfo.GCS,'WGS 84')
            error(['Raster must use the WGS 84 ' ...
                'geographic coordinate system.'])
        end
        if ~strcmpi(sInfo.Datum,'World Geodetic System 1984')
            error('Raster must use the World Geodetic System 1984 datum.')
        end
        if ~strcmpi(sInfo.Ellipsoid,'WGS 84')
            error('Raster must use the WGS 84 ellipsoid.')
        end
        if ~strcmpi(sInfo.PM,'Greenwich')
            error('Raster must use the Greenwich prime meridian.')
        end

    case '.shp'

        % Read .prj file
        try
            strPrj = fileread(strcat([strFileIn(1:end-4) '.prj']));
        catch objExc
            warning(objExc.message)
            [~,strFile,strExt] = fileparts(strFileIn);
            warning(['.prj file not found for ' strFile strExt ...
                '. Assuming the shapefile uses WGS 84 geographic ' ...
                'coordinate system...'])
            return
        end

        % Make sure shapefile uses standard WGS84
        if ~strcmpi(checkShapePrj(strPrj,'GEOGCS["'),'GCS_WGS_1984')
            error(['Shapefile must use the WGS 84 ' ...
                'geographic coordinate system.'])
        end
        if ~strcmpi(checkShapePrj(strPrj,'DATUM["'),'D_WGS_1984')
            error('Shapefile must use the World Geodetic System 1984 datum.')
        end
        if ~any(strcmpi(checkShapePrj(strPrj,'SPHEROID["'), ...
                {'WGS 84','WGS_1984'}))
            error('Shapefile must use the WGS 84 ellipsoid.')
        end
        if ~strcmpi(checkShapePrj(strPrj,'PRIMEM["'),'Greenwich')
            error('Shapefile must use the Greenwich prime meridian.')
        end
        
    otherwise
        error('Invalid file type.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function strVal = checkShapePrj(strPrj,strField)

try
    
    % Find the start index of the field
    iIdxS = strfind(strPrj,strField)+length(strField);

    % Find the end index of the field
    iIdxE = min(iIdxS + strfind(strPrj(iIdxS:end),'"'))-2;
    strVal = strPrj(iIdxS:iIdxE);

catch
    
    % Output empty string
    strVal = '';
    
end
end
end

