function [val paramnames paramnum description]=user2(param,x,opt)

%Beschreibung der Funktion
description='Vielfaches der Fitfunktion bei B=0';

%Anzahl der Parameter:
paramnum=1;

%Namen der Parameter
paramnames = {'scale'};

if (nargin>2) && (strcmp(upper(opt),'INFO'))  val=0; return; end
if (nargin==0) fprintf(['Funktion: ' description]);  paramnames, return; end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

scale=param(1);

p=[0.1262   50.0000    0.1758    0.2848   40.6033    0.1529];
%Aus Fit an therm. Daten 2meV

%p=[0.10764	40.352	0.15531	0.26253	39.961	0.16687]; %4meV-Fit
%p=[0.1621	35.318	0.27199	0.30969	18.816	0.12];    %6meV-Fit
%p=[0.08018	16.209	0.15722	0.26732	25.267	0.2124];  %8meV-Fit


val = scale * (gaussAsym(p(1:3),x,opt) + gaussAsym(p(4:6),x,opt));

