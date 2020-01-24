function mPts = transformUsingSolverVar(mPts,sInput)
% Transform points

% Get variables from input structure
mBnd = sInput.scaleBounds;
vC = sInput.rotationCenter;
vIdx = sInput.variablesIndex;
vVar = sInput.variables;

% Descale the variables
for i = 1:length(vVar)
    vVar(i) = mBnd(1,i) + vVar(i) * (mBnd(2,i) - mBnd(1,i));
end

% Put variables back into correct order
vVarO = [zeros(1,6) ones(1,4)];
vVarO(vIdx) = vVar;

% Compose rotation matrix, translation vector, and scale factors
mR = eye(4); 
mR(1:3,1:3) = makeRotMat(vVarO(1),vVarO(2),vVarO(3),'deg');
vT = [vVarO(4);vVarO(5);vVarO(6)];
vS = [vVarO(7);vVarO(8);vVarO(9)];
dS = vVarO(10);

if strcmpi(sInput.direction,'forward')
    
    % Center the points
    mPts(1,:) = mPts(1,:) - vC(1);
    mPts(2,:) = mPts(2,:) - vC(2);
    mPts(3,:) = mPts(3,:) - vC(3);

    % Apply scale
    mPts(1,:) = mPts(1,:) * vS(1);
    mPts(2,:) = mPts(2,:) * vS(2);
    mPts(3,:) = mPts(3,:) * vS(3);

    % Apply global scale
    mPts(1:3,:) = mPts(1:3,:) * dS;

    % Apply rotation
    mPts = mR * mPts;

    % De-center the points
    mPts(1,:) = mPts(1,:) + vC(1);
    mPts(2,:) = mPts(2,:) + vC(2);
    mPts(3,:) = mPts(3,:) + vC(3);

    % Apply translation
    mPts(1,:) = mPts(1,:) + vT(1);
    mPts(2,:) = mPts(2,:) + vT(2);
    mPts(3,:) = mPts(3,:) + vT(3);

    % Apply polynomial surface correction
    if isfield(sInput,'polySurf')
        if isstruct(sInput.polySurf)
            mPts(3,:) = mPts(3,:) + ...
                polyvaln(sInput.polySurf,[mPts(1,:)' mPts(2,:)'])';
        end
    end

elseif strcmpi(sInput.direction,'inverse')
    
    % Undo polynomial surface correction
    if isfield(sInput,'polySurf')
        if isstruct(sInput.polySurf)
            mPts(3,:) = mPts(3,:) - ...
                polyvaln(sInput.polySurf,[mPts(1,:)' mPts(2,:)'])';
        end
    end
    
    % Undo translation
    mPts(1,:) = mPts(1,:) - vT(1);
    mPts(2,:) = mPts(2,:) - vT(2);
    mPts(3,:) = mPts(3,:) - vT(3);
    
    % Center the points
    mPts(1,:) = mPts(1,:) - vC(1);
    mPts(2,:) = mPts(2,:) - vC(2);
    mPts(3,:) = mPts(3,:) - vC(3);
    
    % Undo rotation
    mPts = mR \ mPts;
    
    % Undo global scale
    mPts(1:3,:) = mPts(1:3,:) / dS;

    % Undo scale
    mPts(1,:) = mPts(1,:) / vS(1);
    mPts(2,:) = mPts(2,:) / vS(2);
    mPts(3,:) = mPts(3,:) / vS(3);
    
    % De-center the points
    mPts(1,:) = mPts(1,:) + vC(1);
    mPts(2,:) = mPts(2,:) + vC(2);
    mPts(3,:) = mPts(3,:) + vC(3);
    
else
    error('Invalid transformation direction parameter.')
end
