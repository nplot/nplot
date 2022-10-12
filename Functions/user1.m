function [val paramnames paramnum description]=user1(param,x,opt)

%Beschreibung der Funktion
description='Function: exp(-ax)*b*x';

%Anzahl der Parameter:
paramnum=2;

%Namen der Parameter
paramnames = {'a', 'b'};

if (nargin>2) && (strcmp(upper(opt),'INFO'))
    val=0;
    return;
end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

val = param(2)*x.*exp(-param(1)*x);
