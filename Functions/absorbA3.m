function [val paramnames paramnum description]=absorb(param,x,opt)

% a3 absorption scan 
% 

description = 'absorb';
paramnames  = {'bgr','mu/ki','x0','A'}; paramnum=4;

if isempty(param),  val=[]; return; end


%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

mu = param(2);
x0 = param(3);

val = param(1)  +  param(4) * exp(-mu ./ cosd(x-x0));
