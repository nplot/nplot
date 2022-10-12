function [qlx,qly,qlz]=QLabA3(a3,gu,gl,a3p,qx,qy,qz,zerovals)

% Q in Lab system
% qx,qy,qz in sample holder system

[stdzeros, GUsign, GLsign] = getoption('stdzeros', 'GUsign', 'GLsign');

% for zero corrections; (if no zeros given, take predefined values)
if nargin<8,   zerovals = stdzeros; end 
[U0, RA30, RGL0, RGU0, RA3P0] = zeromatrix(zerovals);

% Degrees to radians
a3 = pi/180 * a3(:); 
gu = pi/180 * gu(:) * GUsign; 
gl = pi/180 * gl(:) * GLsign;
a3p= pi/180 * a3p(:);



[ok,a3,gu,gl,a3p,qx,qy,qz] = makesamesize(a3,gu,gl,a3p,qx,qy,qz);
if ~ok, return; end


qlx=zeros(size(qx));
qly=qx; qlz=qx;

for i=1:numel(a3)
    RA3 = [cos(a3(i)),  -sin(a3(i)),  0; sin(a3(i)),  cos(a3(i)),  0; 0, 0, 1];
    RGL = [cos(gl(i)), 0, sin(gl(i)); 0, 1, 0; -sin(gl(i)), 0, cos(gl(i))];
    RGU = [1, 0, 0; 0, cos(gu(i)), -sin(gu(i)); 0, sin(gu(i)), cos(gu(i))];
    RA3P = [cos(a3p(i)),  -sin(a3p(i)),  0; sin(a3p(i)),  cos(a3p(i)),  0; 0, 0, 1];
    
    Rs = RA3 * RA30* RGL * RGL0 * RGU * RGU0 * RA3P * RA3P0 * U0';

    for j=1:size(qx,2)
        q = Rs * [-qy(i,j); qx(i,j); qz(i,j)];
        qlx(i,j) =  q(1);
        qly(i,j) =  q(2);
        qlz(i,j) =  q(3);
    end
end
       