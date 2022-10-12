function [val paramnames paramnum description]=lorentzA(param,x,opt)

% lorentzA(param,x)
% Lorentzian peak function (Amplitude version)
% parameters: x0, Amplitude, Fwhm


description='Lorentzian peak';
paramnames = {'x0', 'Amplitude', 'Fwhm'}; paramnum=3;

if isempty(param),  val=[]; return; end

%-------------------------------------------------------
% Definition below
%-------------------------------------------------------

xc    = param(1);
Amplitude = param(2);
FWHM  = param(3);

val = Amplitude * FWHM^2 * 1./(4*(x-xc).^2+FWHM^2);
val = reshape(val,size(x));
end