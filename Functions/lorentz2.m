function [val paramnames paramnum description]=lorentz2(param,x,opt)

%Beschreibung der Funktion
description='Zwei symmetrische Lorentz-Peaks auf konstantem Untergrund';

%Anzahl der Parameter:
paramnum=4;

%Namen der Parameter
paramnames = {'Bgr', 'x0', 'Fläche (ges.)', 'Fwhm'};

if (nargin>2) && (strcmp(upper(opt),'INFO'))
    val=0;
    return;
end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

BG    = param(1);
xc    = param(2);
inten = param(3);
fwhm  = param(4);

val = BG + inten / pi * fwhm *( 1./(4*(x-xc).^2+fwhm^2) + 1./(4*(x+xc).^2+fwhm^2) );
val=reshape(val,size(x));