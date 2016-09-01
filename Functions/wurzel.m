function [val paramnames paramnum description]=wurzel(param,x,opt)

%Beschreibung der Funktion
description='Wurzelfunktion';

%Anzahl der Parameter:
paramnum=2;

%Namen der Parameter
paramnames = {'a', 'x0'};

if (nargin>2) && (strcmp(upper(opt),'INFO'))  val=0; return; end
if (nargin==0) fprintf(['Funktion: ' description]);  paramnames, return; end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

val=param(1)*sqrt(x-param(2));