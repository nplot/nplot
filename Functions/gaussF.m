function [val paramnames paramnum description]=gaussF(param,x,opt)

% gaussA(param,x)
% Gaussian function (Area version)
% parameters: x0, Area, Fwhm

description='Gauss-Peak';
paramnames = {'x0', 'Area', 'Fwhm'};  paramnum=3;

if isempty(param), val=[]; return; end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

xc    = param(1);
inten = param(2);
fwhm  = param(3);


val = inten / sqrt(pi/log(2)) / fwhm * 2 *  exp(-4*log(2)*(x-xc).^2 / fwhm^2);
