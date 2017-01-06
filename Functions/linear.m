function [val, paramnames, paramnum, description]=linear(param,x)
% linear(param,x)
% Linear function f(x) = y0 + slope * x
% Parameters: y0, slope

description = 'linear';
paramnames = {'y0', 'slope'}; paramnum = 2;

if isempty(param),  val=[]; return; end

%-------------------------------------------------------
% Definition below
%-------------------------------------------------------

val = param(1) + x*param(2);
end