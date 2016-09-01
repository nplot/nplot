function [alpha,alphaP]=Scatt_Angle(a4,gfc,ch)

% Calculate true scattering angle and projection of scattering angle onto
% the detector plane
% (replaces a4)
%
% P. Steffens 09/2008

delta = pi/180 * (a4(:)-52.5);
nu    = pi/180 * (15+(2.5*(ch(:)'-1)));
gamma = pi/180 * gfc(:);

cosalpha = cos(delta).*cos(gamma) * cos(nu) - sin(delta) * sin(nu);

alpha = 180/pi * acos(cosalpha);        % (alpha is always positive)

cosalpha31 = -sin(delta);               % alpha for channel 31

cosalphaP31 = cosalpha31 ./  sqrt(cos(gamma).^2  + sin(delta).^2 .* sin(gamma).^2) ;

alphaP31 = 180/pi * acos(cosalphaP31) .* sign(a4(:)+37.5) ;  
                                        % alphaP can be positive or negative
alphaP = alphaP31 * ones(size(nu)) - 2.5 * ones(size(alphaP31)) * (31-ch(:)');


