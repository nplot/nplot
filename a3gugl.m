function [a3,gu,gl,a3p]=a3gugl(a4,gfc,psi,varargin)


% varargin (optional) may contain:
%   opt:      desired solution (give condition, e.g. 'A3P=0', 'A3=25.4' etc.)
%   zerovals: angles zeros
%   config:   spectrometer configuration details (sign_gu etc.)
%
% P. Steffens, 01/2008


delta = pi/180 * (a4(:)-52.5);
gamma = pi/180 *  gfc(:);
psi   = pi/180 *  psi(:);

[ok,delta,gamma,psi] = makesamesize(delta,gamma,psi);
if ~ok, return; end

for m=1:numel(delta)
    
    % Rdelta * Rgamma * Rpsi: (desired rotation matrix)
    DGP=[cos(delta(m))*cos(psi(m))-sin(delta(m))*cos(gamma(m))*sin(psi(m)), -cos(delta(m))*sin(psi(m))-sin(delta(m))*cos(gamma(m))*cos(psi(m)),  sin(delta(m))*sin(gamma(m)) ; ...
         sin(delta(m))*cos(psi(m))+cos(delta(m))*cos(gamma(m))*sin(psi(m)), -sin(delta(m))*sin(psi(m))+cos(delta(m))*cos(gamma(m))*cos(psi(m)), -cos(delta(m))*sin(gamma(m)) ; ...
         sin(gamma(m))*sin(psi(m)),                                          sin(gamma(m))*cos(psi(m)),                                          cos(gamma(m))];
 
    [a3(m),gu(m),gl(m),a3p(m)] = anglesfrommatrix(DGP, varargin{:}); %#ok<AGROW>
end

% a3 = atan2(sin(delta).*cos(psi)+cos(delta).*cos(gamma).*sin(psi), cos(delta).*cos(psi)-sin(delta).*cos(gamma).*sin(psi));
% gu = atan2(sin(gamma).*cos(psi), cos(gamma));
% gl = atan2(-sin(gamma).*sin(psi), sqrt(sin(gamma).^2.*cos(psi).^2+cos(gamma).^2));
% 
% a3 = 180/pi * a3;
% gu = 180/pi * gu * GUsign;
% gl = 180/pi * gl * GLsign;
