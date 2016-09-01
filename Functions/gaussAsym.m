function [val paramnames paramnum description]=gaussAsym(param,x,opt)

% gaussAsym(param,x)
% two Gaussian functions, symmetric around zero (Amplitude version)
% parameters: x0, Amplitude, Fwhm

description='Two symmetric Gauss-Peaks (xc,-xc)';
paramnames = {'xc', 'Amplitude', 'Fwhm'}; paramnum=3;

if isempty(param), val=[]; return; end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

xc   = param(1);
inten= param(2);
fwhm = param(3);

val = inten * ( exp(-4*log(2)*(x-xc).^2 / fwhm^2) +  exp(-4*log(2)*(x+xc).^2 / fwhm^2) );
