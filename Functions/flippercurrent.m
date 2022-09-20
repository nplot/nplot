function [val paramnames paramnum description] = flippercurrent(param,x,opt)

% cosine(param,x,opt)
% Cosine function, f(I) = A * (1-cos (b*(I - Imin))) + BG 
% Parameters: A, b, c

description = 'Cosine';
paramnames = {'A', 'b', 'Imin', 'BG'}; paramnum = length(paramnames);

if isempty(param),  val=[]; return; end

%-------------------------------------------------------
% Function definition
%-------------------------------------------------------

val = param(1) * (1-cos(param(2)*(x - param(3))))/2 + param(4); 
