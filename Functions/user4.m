function [val paramnames paramnum description]=user4(param,x,opt)

%Beschreibung der Funktion
description='user2 verbreitert (Faltung mit Gauss)';

%Anzahl der Parameter:
paramnum=2;

%Namen der Parameter
paramnames = {'scale', 'width'};

if (nargin>2) && (strcmp(upper(opt),'INFO'))  val=0; return; end
if (nargin==0) fprintf(['Funktion: ' description]);  paramnames, return; end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

scale=param(1);
width=param(2);

xgrid=[-4:.1:4];
gauss=exp(-xgrid.^2/2);
gauss=gauss/sum(gauss); %Normierung

val=[];
for xval=x(:)'
    val=[val, sum(user2(scale,xval+width*xgrid,0).*gauss)];
end
val=reshape(val, size(x));
