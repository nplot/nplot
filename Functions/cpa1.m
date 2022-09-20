function [val paramnames paramnum description] = cpa1(param,x,opt)

% Cosine func for CPA alignment

% Enter function name and description, f(x) = A * (1 + cos(2*pi/p * x + x0) + bgr
% Parameters: A, p, x0, bgr

description = 'Cosine';
paramnames = {'Amplitude', 'period', 'x0', 'BGR'}; paramnum = length(paramnames);

if isempty(param),  val=[]; return; end

%-------------------------------------------------------
% Function definition
%-------------------------------------------------------

A=param(1);
p=param(2);
x0=param(3);
bgr=param(4);

val = A * (1 + cos(2*pi/p * x + x0)) + bgr ; % Enter formula
