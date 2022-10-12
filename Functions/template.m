function [val paramnames paramnum description] = template(param,x,opt)

% Template for creation of new function

% Edit and save under new name

% template(param,x,opt)
% Enter function name and description, f(x) = ...
% Parameters: A, lambda, x0

description = '2 Gauss on sloping BG';
paramnames = {'Bgr_Y', 'Bgr_G', 'x0_1', 'Int._1', 'Fwhm_1', 'x0_2', 'Int._2', 'Fwhm_2'}; paramnum = length(paramnames);

if isempty(param),  val=[]; return; end

%-------------------------------------------------------
% Function definition
%-------------------------------------------------------

val = 0; % Enter formula
