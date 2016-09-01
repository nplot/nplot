function [qx,qy,qz]=QSamplePsi(psi,a4,gfc,qlx,qly,qlz)

% Q in cartesian crystal system 
%
% qlx,qly,qlz must be same size; colums correspond to detector channels,
% rows to scan steps (psi,a4,gfc can be single numbers or numel()=number of
% rows in qlx)
%
% P. Steffens 01/2008



psi= psi(:);
a4 = a4(:); 
gfc= gfc(:); %(Factor pi/180 below)
qx=zeros(size(qlx));
qy=qx; qz=qx;

% if (size(psi,1)>1) && (size(psi,1)~=size(qlx,1))
%     disp('Size mismatch between arrays!');
%     return;
% elseif (size(psi,1)==1)
%     psi = psi * ones(size(qlx,1),1);
% end
% 
% if (numel(a4)==1);  a4  = a4  * ones(size(psi)); end
% if (numel(gfc)==1); gfc = gfc * ones(size(psi)); end

[ok,psi,a4,gfc] = makesamesize(psi,a4,gfc,qlx(:,1),qly(:,1),qlz(:,1));
if ~ok, return; end

delta =-pi/180*(90-37.5-a4);
gamma = pi/180*gfc;
psi   = pi/180*psi;
    
for i=1:numel(psi)
    
    RG = [1, 0, 0;  0, cos(gamma(i)), -sin(gamma(i));  0, sin(gamma(i)), cos(gamma(i))];
    RD = [cos(delta(i)), -sin(delta(i)), 0; sin(delta(i)), cos(delta(i)), 0; 0, 0, 1];
    RPsi0 = [cos(psi(i)), -sin(psi(i)), 0; sin(psi(i)), cos(psi(i)), 0; 0, 0, 1];

    
    %Rs = RD * RG * RPis0;
    Rs_inv = RPsi0' * RG' * RD';

    for j=1:size(qlx,2)
        q = Rs_inv * [qlx(i,j); qly(i,j); qlz(i,j)];
        qx(i,j) = q(2);
        qy(i,j) =-q(1);
        qz(i,j) = q(3);
    end
end
        
