function [val paramnames paramnum description] = template(param,x,opt)

% Template for creation of new function

% Edit and save under new name

% template(param,x,opt)
% Enter function name and description, f(x) = ...
% Parameters: A, lambda, x0

description = 'no name';
paramnames = {'P1', 'P2', '...'}; paramnum = length(paramnames);

if isempty(param),  val=[]; return; end

%-------------------------------------------------------
% Function definition
%-------------------------------------------------------

val = 0; % Enter formula
