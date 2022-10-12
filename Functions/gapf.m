function [val paramnames paramnum description] = gapf(param,x,opt)

% gapf(param,x)
% guide-to-the-eye function for SC gap
% parameters: a,b,c,d

description  = 'Gap function';
paramnames   = {'a', 'b', 'c', 'd'};  paramnum = length(paramnames);

if isempty(param) || ((nargin>2) && (strcmpi(opt,'INFO'))),  val=[]; return; end

%-------------------------------------------------------
% function definition
%-------------------------------------------------------

a  = param(1);
b  = param(2);
c  = param(3);
d  = param(4);


val = a * (1+tanh(sqrt(b*x)-c)) .* exp(-x*d);
