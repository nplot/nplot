function [qx,qy]=inplaneQ(xang,yang,ki,kf,Qvert)

% Calculate the in-plane components of the scattering vector
%
% Qvert (vertical component of Q), ki, kf all fixed single values or same
% size as xyang
% xang, yang as defined elsewhere (xang projection of Scattering angle on
% the plane, yang=psi+xang(ch31)-90)
%
% P. Steffens, 01/2008

% xx_par denotes projections

alphaP= pi/180 * xang;
beta  = pi/180 * yang;

ki_par = sqrt(ki.^2-Qvert.^2);

Q_par_x = kf.*cos(alphaP) - ki_par;
Q_par_y = kf.*sin(alphaP);

qx =  cos(beta) .* Q_par_x + sin(beta) .* Q_par_y;
qy = -sin(beta) .* Q_par_x + cos(beta) .* Q_par_y;

QSign = getoption('QSign');
if QSign==-1
    qx = -qx; qy = -qy;
end