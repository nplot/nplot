function [val paramnames paramnum description]=gaussFsym(param,x,opt)

% gaussAsym(param,x)
% two Gaussian functions, symmetric around zero (Area version)
% parameters: x0, Area (one peak), Fwhm


description = 'Two symmetric Gaussian peaks (x0,-x0)';
paramnames = {'x0', 'Area(each)', 'Fwhm'}; paramnum=3;

if isempty(param), val=[]; return; end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

xc=param(1);
inten=param(2);
fwhm=param(3);


val = inten / sqrt(pi/log(2)) / fwhm * 2 * ( exp(-4*log(2)*(x-xc).^2 / fwhm^2) + exp(-4*log(2)*(x+xc).^2 / fwhm^2) );
