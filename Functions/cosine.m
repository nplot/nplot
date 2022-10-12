function [val paramnames paramnum description] = cosine(param,x,opt)

% cosine(param,x,opt)
% Cosine function, f(x) = A * cos (b*x + c) (in radians)
% Parameters: A, b, c

description = 'Cosine';
paramnames = {'A', 'b', 'c'}; paramnum = length(paramnames);

if isempty(param),  val=[]; return; end

%-------------------------------------------------------
% Function definition
%-------------------------------------------------------

val = param(1) * cos(param(2)*x + param(3)); 
