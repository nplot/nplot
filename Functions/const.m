function [val paramnames paramnum description]=const(param,x,opt)

% const (param, x)
% Contant value f(x) = const  (=param(1))

description = 'constant';
paramnames  = {'const'}; paramnum=1;

if isempty(param),  val=[]; return; end


%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

val = param(1) * ones(size(x));
