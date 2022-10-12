function [val paramnames paramnum description]=lorentzSpec(param,x,opt)

%Beschreibung der Funktion
description='Lorentz-Fkt. symm. (Spektrum)';

%Anzahl der Parameter:
paramnum=3;

%Namen der Parameter
paramnames = {'x0', 'Amplitude', 'Fwhm'};

if (nargin>2) && (strcmp(upper(opt),'INFO'))  val=0; return; end
if (nargin==0) fprintf(['Funktion: ' description]);  paramnames, return; end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

xc    = param(1);
inten = param(2);
fwhm  = param(3);

val = inten * fwhm^2 * ( 1./(4*(x-xc).^2+fwhm^2) + 1./(4*(x+xc).^2+fwhm^2)  ) ;
val=reshape(val,size(x));