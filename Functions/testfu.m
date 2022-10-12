function [val paramnames paramnum description] = testfu(x,param)

% gaussA(param,x)
% Gaussian function
% parameters: x0, Amplitude, Fwhm

description  = 'Gauss-Peak';
paramnames   = {'x0', 'Amplitude', 'Fwhm'};  paramnum = length(paramnames);

if isempty(param) || ((nargin>2) && (strcmpi(opt,'INFO'))),  val=[]; return; end

%-------------------------------------------------------
% function definition
%-------------------------------------------------------

% xc    = param(2);
% inten = param(1);
% fwhm  = param(3);
% 
% 
% val = inten * exp(-4*log(2)*(x-xc).^2 / fwhm^2) + param(4);

val = param(1) * ones(size(x));