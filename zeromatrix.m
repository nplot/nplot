function [U0, RA30, RGL0, RGU0, RA3P0] = zeromatrix(zerovals, config)

% Matrices for zeros corrections...

% if second parameter given, save call to getoption
if nargin < 2
    [za3_0, zgu_0, zgl_0, GUsign, GLsign] = getoption('za3_0', 'zgu_0', 'zgl_0', 'GUsign', 'GLsign');
else
    za3_0 = config.za3_0; zgu_0 = config.zgu_0; zgl_0 = config.zgl_0;
    GUsign = config.GUsign; GLsign = config.GLsign;
end

gu0 = - pi/180 * (zerovals.gu-zgu_0) * GUsign; 
gl0 = - pi/180 * (zerovals.gl-zgl_0) * GLsign; 
a30 = - pi/180 * (zerovals.a3-za3_0); 
a3p0= - pi/180 * zerovals.a3p;
% "-" signs here because a shift dx in zgu corresponds to a physical
% desplacement -dx for same gu-value


RA30 = [cos(a30),  -sin(a30),  0; sin(a30),  cos(a30),  0; 0, 0, 1];
RGL0 = [cos(gl0), 0, sin(gl0); 0, 1, 0; -sin(gl0), 0, cos(gl0)];
RGU0 = [1, 0, 0; 0, cos(gu0), -sin(gu0); 0, sin(gu0), cos(gu0)];
RA3P0= [cos(a3p0),  -sin(a3p0),  0; sin(a3p0),  cos(a3p0),  0; 0, 0, 1];
U0   = RA30 * RGL0 * RGU0 * RA3P0;