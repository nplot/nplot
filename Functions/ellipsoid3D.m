function [val paramnames paramnum description]=ellipsoid3D(param,x,opt)

%Beschreibung der Funktion
description='Ellipsoid im 3D Raum';

%Anzahl der Parameter:
paramnum=7;

%Namen der Parameter
paramnames = {'a11', 'a12', 'a13', 'a22', 'a23', 'a33', 'Scale'};

if (nargin>2) && (strcmp(upper(opt),'INFO'))
    val=0;
    return;
end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

A= [param(1), param(2), param(3); ...
    param(2), param(4), param(5); ...
    param(3), param(5), param(6)];

for i=1:size(x,1)
    val(i)  = param(7) * exp(- x(i,:) * A * x(i,:)');
end

val = val(:);
