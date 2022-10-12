function [val paramnames paramnum description]=absorb(param,x,opt)

% const (param, x)
% Contant value f(x) = const  (=param(1))

description = 'absorb';
paramnames  = {'bgr','mu','x0','A1','A2','ki1','ki2'}; paramnum=7;

if isempty(param),  val=[]; return; end


%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

mu = param(2);
x0 = param(3);



a3 = x(:,1);
ki(x(:,2)==1) = param(6);
ki(x(:,2)==2) = param(7);

A(x(:,2)==1) = param(4);
A(x(:,2)==2) = param(5);

val = param(1)  +  A' .* exp(-mu./ki' ./ cosd(a3-x0));
