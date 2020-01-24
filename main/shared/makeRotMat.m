function mR = makeRotMat(dRx,dRy,dRz,strOpt)
% Make rotation matrix from 3 euler angles

% Convert to radians if necessary  
if strcmp(strOpt,'deg') == 1
    dRx = dRx * pi/180;
    dRy = dRy * pi/180;
    dRz = dRz * pi/180;
else
    error('specify "deg" or "rad" as 4th input.')
end

% Rotation around x axis
mRx = [1 0         0
       0 cos(dRx) -sin(dRx)
       0 sin(dRx)  cos(dRx)];

% Rotation around y axis
mRy = [cos(dRy) 0 sin(dRy)
       0        1 0
      -sin(dRy) 0 cos(dRy)];

% Rotation around z axis
mRz = [cos(dRz) -sin(dRz) 0
       sin(dRz)  cos(dRz) 0
       0         0        1];

% Multiply rotation matrices
mR = mRx * mRy * mRz;