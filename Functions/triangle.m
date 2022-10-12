function [val paramnames paramnum description] = triangle(param,x,opt)

% Triangle peak function

% triangle(param,x,opt)
% Enter function name and description, f(x) = triang(x)
% Parameters: pos, height, width

description = 'triangle peak';
paramnames = {'x0', 'h', 'fwhm'}; paramnum = length(paramnames);

if isempty(param),  val=[]; return; end

%-------------------------------------------------------
% Function definition
%-------------------------------------------------------

x0 = param(1);
h  = param(2);
fwhm=param(3);

% val(x <= x0-fwhm | x >= x0+fwhm) = 0; 
val = zeros(size(x));
val(x > x0-fwhm & x <= x0) = h/fwhm * (-x0+fwhm + x(x > x0-fwhm & x <= x0)); 
val(x > x0 & x < x0+fwhm)  = h/fwhm *  (x0+fwhm - x(x > x0 & x < x0+fwhm));

