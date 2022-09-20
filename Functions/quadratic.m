function [val paramnames paramnum description] = quadratic(param,x,opt)

% Template for creation of new function

% Edit and save under new name

% template(param,x,opt)
% Enter function name and description, f(x) = ...
% Parameters: A, lambda, x0

description = 'a*(x-x0)^2+y0';
paramnames = {'x0','a','y0'}; paramnum = length(paramnames);

if isempty(param),  val=[]; return; end

%-------------------------------------------------------
% Function definition
%-------------------------------------------------------

x0 = param(1);
a  = param(2);
b  = param(3);

val = a*(x-x0).^2 + b; % Enter formula
