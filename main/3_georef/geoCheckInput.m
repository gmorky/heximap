function [] = geoCheckInput(strRef,strShpPath,hW)

% Update waitbar
try
waitbar(1/8,hW,'checking user input for errors...')
catch
end

% Check reference DEM
checkInput(strRef);

% Check shapefiles
checkShpPath(strShpPath);
