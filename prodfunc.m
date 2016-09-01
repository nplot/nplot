function [val paramnames paramnum description]=prodfunc(param,x,varargin)

%Hilfsfunktion für Produktbildung; wird nicht direkt aufgerufen 
%("fmult" benutzen!)


paramnum=0;
description='(';


%Beschreibung der Funktion und Parameter
for c=1:(nargin-2)
    if c>1, description = [description ' * ']; end
    func = varargin{c};
    [erg names pnum(c) des] = func([],[]);
    description = [description des];
    for k = 1:pnum(c)
        paramnames{paramnum+k} = ['(' num2str(c) ')' names{k}];
    end
    paramnum = paramnum + pnum(c);
end
description = [description ')'];

if isempty(param), val=[]; return; end

%-------------------------------------------------------
% Ab hier Definition der Funktion
%-------------------------------------------------------

val=1;
p=0;
for c = 1:(nargin-2)
    func = varargin{c};
    val = val .* func(param(p+1:p+pnum(c)),x);
    p = p + pnum(c);
end
