function [qlx,qly,qlz]=QLab(a4,gfc,ch,ki,kf)

% Q-vector in Laboratory coordinates.
% columns of output correspond to detector channels,
% rows to scan steps (a4,gfc,ki,kf same length or single numbers)
%
% P. Steffens 01/2008


a4=a4(:); gfc=gfc(:); ki=ki(:); kf=kf(:);
ch=ch(:)';

delta =-pi/180*(90-37.5-a4);
gamma = pi/180*gfc;
nu    = pi/180*(15+(2.5*(ch-1)));

qlx = - (kf .* cos(delta)) * sin(nu) - (kf .* sin(delta) .* cos(gamma)) * cos(nu);
qly = - (kf .* sin(delta)) * sin(nu) + (kf .* cos(delta) .* cos(gamma)) * cos(nu)  - ki * ones(size(nu));
qlz =   (kf .* sin(gamma)) * cos(nu);

% DGN=[[cos(delta)*cos(nu)-sin(delta)*cos(gamma)*sin(nu), -cos(delta)*sin(nu)-sin(delta)*cos(gamma)*cos(nu), sin(delta)*sin(gamma)]; ...
%      [sin(delta)*cos(nu)+cos(delta)*cos(gamma)*sin(nu), -sin(delta)*sin(nu)+cos(delta)*cos(gamma)*cos(nu), -cos(delta)*sin(gamma)]; ...
%      [sin(gamma)*sin(nu), sin(gamma)*cos(nu), cos(gamma)]];
 
QSign = getoption('QSign');
% if Q=ki-kf desired, change signs
if QSign==-1
    qlx = -qlx;  qly = -qly;  qlz = -qlz;
end

