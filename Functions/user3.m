function [val paramnames paramnum description]=user3(param,x,opt)

%Beschreibung der Funktion
description='Zwei um Null symmetrische Gauss-Peaks mit versch. FWHM';

%Anzahl der Parameter:
paramnum=4;

%Namen der Parameter
paramnames = {'x0', 'Fläche(jew.)', 'Fwhm1', 'Fwhm2'};

if (nargin>2) && (strcmp(upper(opt),'INFO'))  val=0; return; end
if (nargin==0) fprintf(['Funktion: ' description]);  paramnames, return; end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

xc=param(1);
inten=param(2);
fwhm1=param(3);
fwhm2=param(4);


val = inten / sqrt(pi/log(2))  * 2 * ( exp(-4*log(2)*(x-xc).^2 / fwhm1^2)/fwhm1 + exp(-4*log(2)*(x+xc).^2 / fwhm2^2)/fwhm2 );
