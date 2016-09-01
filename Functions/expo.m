function [val paramnames paramnum description] = expo(param,x,opt)

% expo(param,x,opt)
% Exponential function f(x) = A*exp((x-x0)*lambda)
% Parameters: A, lambda, x0

description='Exponentialfunktion[A*exp((x-x0)*l)]';
paramnames = {'A', 'lambda', 'x0'}; paramnum = length(paramnames);

if isempty(param),  val=[]; return; end

%-------------------------------------------------------
% Function definition
%-------------------------------------------------------

val = param(1) * exp((x-param(3))*param(2));
