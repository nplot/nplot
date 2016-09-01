function [val paramnames paramnum description]=gauss(param,x,opt)

%Beschreibung der Funktion
description='Hyperbel a/(x^e-b)';

%Anzahl der Parameter:
paramnum=3;

%Namen der Parameter
paramnames = {'a', 'b', 'e'};

if (nargin>2) && (strcmp(upper(opt),'INFO'))
    val=0;
    return;
end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

a=param(1);
b=param(2);
e=param(3);

val = a ./ (x.^e-b);
