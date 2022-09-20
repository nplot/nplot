function [val paramnames paramnum description] = cosine_2(param,x,opt)

% cosine(param,x,opt)
% Cosine function, f(x) = A * cos (b*(x - x0)) + C (in radians)
% Parameters: A, b, c

description = 'Cosine';
paramnames = {'A', 'b', 'x0', 'C'}; paramnum = length(paramnames);

if isempty(param),  val=[]; return; end

%-------------------------------------------------------
% Function definition
%-------------------------------------------------------

val = param(1) * cos(param(2)*(x - param(3))) + param(4); 
