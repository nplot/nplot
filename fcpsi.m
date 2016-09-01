function p=fcpsi(a3,gu,gl,a3p,a4,gfc,zerovals)

% Calculate psi from the other spectrometer angles
%
% P. Steffens 01/2008



if nargin < 7, zerovals = getoption('stdzeros'); end

% a3 = pi/180*a3(:);
% gu = pi/180*gu(:) * GUsign;
% gl = pi/180*gl(:) * GLsign;
% a3p= pi/180*a3p(:);
gamma = pi/180*gfc(:);
delta = pi/180*(a4(:)-52.5);

psiref_1 = -sin(delta) .* cos(gamma);
psiref_2 =  cos(delta) .* cos(gamma);
psiref_3 =      sin(gamma);

% a_1 = cos(a3).*sin(gl).*sin(gu) - sin(a3).*cos(gu);
% a_2 = sin(a3).*sin(gl).*sin(gu) + cos(a3).*cos(gu);
% a_3 =    cos(gl).*sin(gu);

[a_1,a_2,a_3] = QLabA3(a3,gu,gl,a3p,1,0,0,zerovals);

cospsi = psiref_1.*a_1 + psiref_2.*a_2 + psiref_3.*a_3;

p = 180/pi* real(acos(cospsi)) .* sign(.1+sign(psiref_1.*a_2 - psiref_2.*a_1));