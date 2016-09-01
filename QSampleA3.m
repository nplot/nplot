function [qx,qy,qz]=QSampleA3(a3,gu,gl,a3p,qlx,qly,qlz,zerovals)

% Q in cartesian crystal system
% qlx,qly,qlz in Lab system
%
% P. Steffens 02/2008

%options;
[stdzeros, GUsign, GLsign] = getoption('stdzeros', 'GUsign', 'GLsign');

% for zero corrections; (if no zeros given, take predefined values)
if nargin<8,   zerovals = stdzeros; end 
[U0, RA30, RGL0, RGU0, RA3P0] = zeromatrix(zerovals);

% Degrees to radians
a3 = pi/180 * a3(:); 
gu = pi/180 * gu(:) * GUsign; 
gl = pi/180 * gl(:) * GLsign;
a3p= pi/180 * a3p(:);

qx=zeros(size(qlx));
qy=qx; qz=qx;

[ok,a3,gu,gl,a3p] = makesamesize(a3,gu,gl,a3p);
if ~ok, return; end
    
for i=1:numel(a3)
    RA3 = [cos(a3(i)),  -sin(a3(i)),  0; sin(a3(i)),  cos(a3(i)),  0; 0, 0, 1];
    RGL = [cos(gl(i)), 0, sin(gl(i)); 0, 1, 0; -sin(gl(i)), 0  cos(gl(i))];
    RGU = [1, 0, 0; 0, cos(gu(i)), -sin(gu(i)); 0, sin(gu(i)), cos(gu(i))];
    RA3P = [cos(a3p(i)),  -sin(a3p(i)),  0; sin(a3p(i)),  cos(a3p(i)),  0; 0, 0, 1];
    
    % Rs = RA3 * RA30* RGL * RGL0 * RGU * RGU0 * RA3P * RA3P0 * U0';
    Rs_inv =  U0 * RA3P0' * RA3P' * RGU0' * RGU' * RGL0' * RGL' * RA30' * RA3';

    for j=1:size(qlx,2)
        q = Rs_inv * [qlx(i,j); qly(i,j); qlz(i,j)];
        qx(i,j) =  q(2);
        qy(i,j) = -q(1);
        qz(i,j) =  q(3);
    end
end
        
