function [h,k,l]=calcHKL(qx,qy,qz,UB)

% Calculate hkl indices for scattering vector (qx,qy,qz) in crystal
% cartesian coordinates.

[ok,qx,qy,qz] = makesamesize(qx,qy,qz);

h = zeros(size(qx)); 
k = zeros(size(qx)); 
l = zeros(size(qx)); 


% qx=qx(:);
% qy=qy(:);
% qz=qz(:);

for i=1:numel(qx)
    hkl = UB \ [qx(i); qy(i); qz(i)];  
    % UB \ x means UB^(-1)*x 
    h(i)=hkl(1);
    k(i)=hkl(2);
    l(i)=hkl(3);
end

% h = reshape(h,size(qx));
% k = reshape(k,size(qy));
% l = reshape(l,size(qz));
