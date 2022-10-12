function rot = rotationtoviewingplane(ax)

% Get the Camera view parameters of the axes "ax" to construct a right
% handed normal system with the z-axis pointing towards the camera
% rot is the rotation matrix to be applied to axis coordinates

% P. Steffens 08/2014

camPos    = get(ax, 'CameraPosition');  % camera position
camTgt    = get(ax, 'CameraTarget');    % where the camera is pointing to
camUpVect = get(ax, 'CameraUpVector');  % camera 'up' vector
camDir = camPos - camTgt; % camera direction

% build an orthonormal frame based on the viewing direction and the 
% up vector (the "view frame")
upAxis= camUpVect/norm(camUpVect); 
zAxis = camDir/norm(camDir);    
xAxis = cross(upAxis, zAxis);
yAxis = cross(zAxis, xAxis);

rot = [xAxis; yAxis; zAxis]; % Rotation matrix