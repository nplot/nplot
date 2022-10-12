function [a1,a2,b] = fitplane(x1,x2,y)

x1=x1(:);
x2=x2(:);

if nargin == 2
    mat = [sum(x1.^2), sum(x1); sum(x1), numel(x1)];
    erg = mat^(-1) * [sum(x1.*x2); sum(x2)];
    a1 = erg(1); a2=erg(2);
    
else


y=y(:);

mat = [ sum(x1.^2),   sum(x1.*x2), sum(x1); ...
        sum(x1.*x2),  sum(x2.^2),  sum(x2); ...
        sum(x1),      sum(x2),     numel(x1)];
    
erg = inv(mat) * [ sum(x1.*y); sum(x2.*y); sum(y)];

a1 = erg(1);
a2 = erg(2);
b  = erg(3);

end         