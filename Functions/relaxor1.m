function [val paramnames paramnum description]=relaxor1(param,x,opt)

%Beschreibung der Funktion
description='"Single Relaxor" mit lin. BG';

%Anzahl der Parameter:
paramnum=5;

%Namen der Parameter
paramnames = {'chi0', 'Gamma', 'w0', 'BG_0', 'BG_lin'};

if (nargin>2) && (strcmp(upper(opt),'INFO'))  val=0; return; end
if (nargin==0) fprintf(['Funktion: ' description]);  paramnames, return; end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

chi0=param(1);
Gamma=param(2);
w0=param(3);

val = chi0 * x * Gamma ./ ((x-w0).^2+Gamma^2) + param(4) + param(5)*x; ;
