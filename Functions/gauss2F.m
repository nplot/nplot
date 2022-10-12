function [val paramnames paramnum description]=gauss2F(param,x,opt)

%Beschreibung der Funktion
description='Gauss-Funktion mit konstantem Untergrund';

%Anzahl der Parameter:
paramnum=7;

%Namen der Parameter
paramnames = {'Bgr', 'x0_1', 'Fläche_1', 'Fwhm_1', 'x0_2', 'Fläche_2', 'Fwhm_2'};

if (nargin>2) && (strcmp(upper(opt),'INFO'))
    val=0;
    return;
end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

BG=param(1);
xc1=param(2);
inten1=param(3);
fwhm1=param(4);
xc2=param(5);
inten2=param(6);
fwhm2=param(7);

val = BG + inten1 / sqrt(pi/log(2)) / fwhm1 * 2 *  exp(-4*log(2)*(x-xc1).^2 / fwhm1^2) + inten2 / sqrt(pi/log(2)) / fwhm2 * 2 *  exp(-4*log(2)*(x-xc2).^2 / fwhm2^2) ;
