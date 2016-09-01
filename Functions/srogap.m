function [val, paramnames, paramnum, description]=srogap(param,x)
% srogap(param,x)

description = 'gapfunc';
paramnames = {'chi0', 'Gamma','Egap','width'}; paramnum = 4;

if isempty(param),  val=[]; return; end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

val = param(1)* x * param(2) ./ (param(2)^2+x.^2) .* (tanh((x-param(3))*pi/param(4))+1)/2  ;
