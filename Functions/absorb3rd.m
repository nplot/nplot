function [val paramnames paramnum description]=absorb3rd(param,x,opt)

% const (param, x)
% Contant value f(x) = const  (=param(1))

description = 'absorb';
paramnames  = {'bgr','mu','x0','A1','A2','A3','ki_nom'}; paramnum=7;

if isempty(param),  val=[]; return; end


%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

mu = param(2);
x0 = param(3);

ki_nom = param(7);
A1 = param(4);
A2 = param(5);
A3 = param(6);


a3 = x;

val = param(1)  +  A1/ki_nom * exp(-mu/ki_nom ./ cosd(a3-x0)) + A2/2/ki_nom * exp(-mu/ki_nom/2 ./ cosd(a3-x0)) + A3/3/ki_nom * exp(-mu/ki_nom/3 ./ cosd(a3-x0))  ;
