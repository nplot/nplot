function [val paramnames paramnum description] = gaussAcut(param,x,opt)

% gaussA(param,x)
% Gaussian function (Amplitude version)
% parameters: x0, Amplitude, Fwhm, cutparam

description  = 'Gaussian peak';
paramnames   = {'x0', 'Amplitude', 'Fwhm', 'Cutoff'};  paramnum = length(paramnames);

if isempty(param) || ((nargin>2) && (strcmpi(opt,'INFO'))),  val=[]; return; end

%-------------------------------------------------------
% function definition
%-------------------------------------------------------

xc    = param(1);
inten = param(2);
fwhm  = param(3);
cutoff = param(4);


val = cutoff * (atan(inten * exp(-4*log(2)*(x-xc).^2 / fwhm^2) / cutoff));
