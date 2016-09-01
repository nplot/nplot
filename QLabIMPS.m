function [qlx,qly,qlz]=QLabIMPS(a4,a5,ki,kf,roi,lthree)

% Q-vector in Laboratory coordinates.

% P. Steffens 03/2011


a4=a4(:); a5=a5(:); ki=ki(:); kf=kf(:); roi = roi(:);

 
[QSign, za] = getoption('QSign','impsanadist');

delta    = atand( (roi-5) * za .* sind(a5) ./ ( (roi-5)* za .* cosd(a5) - lthree*100 ) );

twotheta = a4 + delta;


qlx = - kf .* sind(twotheta); 
qly =   kf .* cosd(twotheta)  - ki ;
qlz =   zeros(size(kf));

% if Q=ki-kf desired, change signs
if QSign==-1
    qlx = -qlx;  qly = -qly;  qlz = -qlz;
end

