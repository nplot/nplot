function b=bose(E,T,opt)

if (nargin==0)
    fprintf('Bosefaktor. Syntax: bose(E, T [,´THz´])\n');
    return
end
b=0;
if (nargin>2) && (strcmp(upper(opt),'THZ'))
        E=E*4.135;
end
b=1./(1-exp(-E/0.08617./T));
