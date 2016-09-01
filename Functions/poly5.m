function [val paramnames paramnum description]=poly5(param,x,opt)

%Beschreibung der Funktion
description='Polynom 5. Grades';

%Anzahl der Parameter:
paramnum=6;

%Namen der Parameter
paramnames = {'a0', 'a1', 'a2', 'a3', 'a4', 'a5'};

if isempty(param) | ((nargin>2) & (strcmpi(opt,'INFO'))),  val=[]; return; end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

val=0;
for k=1:1:paramnum
    val = val + param(k) * x.^(k-1);
end
