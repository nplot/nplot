function [a4new, psioffset] = a4step(a4,gfc,apstep)

% [a4new, psioffset] = a4step(a4,gfc,apstep)
% Calculate the new a4-angle starting from the present one and the desired
% change in the projection of the scattering angle onto the plane (alpha').
% Also calulate the offset of the new psi-scan to be well aligned with the previous one.  
%
% a4,gfc: current values
%
% P. Steffens, 03/2009

delta = pi/180 * (a4-52.5);
gamma = pi/180 * gfc;
x = pi/180*apstep;

sindeltanew = cos(x)*sin(delta) + sin(x)*cos(delta)*cos(gamma);

a4new = 180/pi * asin(sindeltanew) + 52.5;

cosa31p_new = -sindeltanew / sqrt(1 - cos(delta)*sin(gamma));
cosa31p_old = -sin(delta) / sqrt(1 - cos(delta)*sin(gamma));

psioffset = 180/pi * (acos(cosa31p_old) - acos(cosa31p_new));

gamma_new = 180/pi * asin(cos(delta)*sin(gamma)/sqrt(1-sindeltanew^2));

fprintf([' To change the detector angle by %2.2f degrees, the new a4 value is %2.2f .\n' ...
    ' The new gfc will be %2.2f .\n' ...
    ' The next psi-scan should be shifted by %2.2f to be well aligned with the previous ones.\n'], apstep, a4new, gamma_new, psioffset);
