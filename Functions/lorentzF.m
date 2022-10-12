function [val paramnames paramnum description]=lorentzF(param,x,opt)

% lorentzF(param,x)
% Lorentzian peak function (Area version)
% parameters: x0, Area, Fwhm

description='Lorentzian peak';
paramnames = {'x0', 'Fläche', 'Fwhm'}; paramnum=3;

if isempty(param),  val=[]; return; end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

xc    = param(1);
inten = param(2);
fwhm  = param(3);

val = inten *2/pi * fwhm * 1./(4*(x-xc).^2+fwhm^2);
