function [val paramnames paramnum description]=gaussbg(param,x,opt)

%Beschreibung der Funktion
description='Gauss-Funktion mit konstantem Untergrund';

%Anzahl der Parameter:
paramnum=4;

%Namen der Parameter
paramnames = {'Bgr', 'x0', 'Amplitude', 'Fwhm'};

if (nargin>2) && (strcmp(upper(opt),'INFO'))
    val=0;
    return;
end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

BG=param(1);
xc=param(2);
inten=param(3);
fwhm=param(4);


val = BG + inten * exp(-4*log(2)*(x-xc).^2 / fwhm^2);
