function [val paramnames paramnum description] = pvoigtcut(param,x,opt)

% pvoigt(param,x)
% Pseudo-Voigt function 
% parameters: x0, Amplitude, width, eta

description  = 'Pseudo-Voigt';
paramnames   = {'x0', 'Amplitude', 'Width', 'Eta', 'Cutoff'};  paramnum = length(paramnames);

if isempty(param) || ((nargin>2) && (strcmpi(opt,'INFO'))),  val=[]; return; end

%-------------------------------------------------------
% function definition
%-------------------------------------------------------

xc    = param(1);
inten = param(2);
fwhm  = param(3);
eta   = param(4);
cutoff= param(5);

val = cutoff * atan(1/cutoff* inten * ( eta * lorentzA([xc,1,fwhm],x) + (1-eta) * gaussA([xc,1,fwhm],x) ));

end

