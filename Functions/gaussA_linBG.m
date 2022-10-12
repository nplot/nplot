function [val paramnames paramnum description]=gaussA_linBG(param,x,opt)

%Beschreibung der Funktion
description='Gauss-Funktion mit linearem Untergrund';

%Anzahl der Parameter:
paramnum=5;

%Namen der Parameter
paramnames = {'Bgr_const', 'Bgr_lin', 'x0', 'Amplitude', 'Fwhm'};

if (nargin>2) && (strcmp(upper(opt),'INFO'))
    val=0;
    return;
end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

BG=param(1);
slope=param(2);
xc=param(3);
inten=param(4);
fwhm=param(5);


val = BG + slope*x + inten * exp(-4*log(2)*(x-xc).^2 / fwhm^2);
