function [qx,qy,qz]=calcQS(h,k,l,UB)

% Calculate Q in sample holder system from HKL's

h=h(:);
k=k(:);
l=l(:);

for i=1:numel(h)
    qq = UB * [h(i); k(i); l(i)];  
    qx(i)=qq(1);
    qy(i)=qq(2);
    qz(i)=qq(3);
end

qx = reshape(qx,size(h));
qy = reshape(qy,size(k));
qz = reshape(qz,size(l));

