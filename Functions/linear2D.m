function [val paramnames paramnum description]=linear2D(param,x)
% linear(param,x)
% Linear function f(x) = y0 + a1 * x(:,1) + a2 * x(:,2)
% Parameters: y0, a1, a2


description = 'linear';
paramnames = {'y0', 'a1', 'a2'}; paramnum = 3;

if isempty(param),  val=[]; return; end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

val = param(1) + x(:,1)*param(2) + x(:,2)*param(3);
