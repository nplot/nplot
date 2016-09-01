function [val paramnames paramnum description] = pvoigt(param,x,opt)

% pvoigt(param,x)
% Pseudo-Voigt function 
% parameters: x0, Amplitude, width, eta

description  = 'Pseudo-Voigt';
paramnames   = {'x0', 'Amplitude', 'Width', 'Eta'};  paramnum = length(paramnames);

if isempty(param) || ((nargin>2) && (strcmpi(opt,'INFO'))),  val=[]; return; end

%-------------------------------------------------------
% function definition
%-------------------------------------------------------

xc    = param(1);
inten = param(2);
fwhm  = param(3);
eta   = param(4);

val = inten * ( eta * lorentzA([xc,1,fwhm],x) + (1-eta) * gaussA([xc,1,fwhm],x) );
