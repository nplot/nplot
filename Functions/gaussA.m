function [val, paramnames, paramnum, description] = gaussA(param,x,opt)

% gaussA(param,x)
% Gaussian function (Amplitude version)
% parameters: x0, Amplitude, Fwhm

description  = 'Gaussian peak';
paramnames   = {'x0', 'Amplitude', 'Fwhm'};  paramnum = length(paramnames);

if isempty(param) || ((nargin>2) && (strcmpi(opt,'INFO'))),  val=[]; return; end

%-------------------------------------------------------
% function definition
%-------------------------------------------------------

xc        = param(1);
Amplitude = param(2);
FWHM      = param(3);


val = Amplitude * exp(-4*log(2)*(x-xc).^2 / FWHM^2);
end