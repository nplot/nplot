function [val paramnames paramnum description]=lorentzA(param,x,opt)

% lorentzA(param,x)
% Lorentzian peak function (Amplitude version)
% parameters: x0, Amplitude, Fwhm


description='Lorentzian peak';
paramnames = {'x0', 'Amplitude', 'Fwhm'}; paramnum=3;

if isempty(param),  val=[]; return; end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

xc    = param(1);
inten = param(2);
fwhm  = param(3);

val = inten * fwhm^2 * 1./(4*(x-xc).^2+fwhm^2);
val = reshape(val,size(x));