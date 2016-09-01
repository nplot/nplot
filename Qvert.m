function Qv = Qvert(scan, nomean)

% for FLATCONE:
% Determine vertical component of momentum transfer
% If nomean, no averaging done
%
% P. Steffens 07/2008

[QSign, maxdeviate] = getoption('QSign','maxdeviate');

delta =-pi/180*(90-37.5-getvar(scan,'a4'));
gamma = pi/180*getvar(scan,'gfc');
ki = getvar(scan,'ki');

Qv = cos(delta).*sin(gamma).*ki;

if nargin<2 || ~nomean  % Give average unless nomean is set
    if (max(Qv)-min(Qv)) < maxdeviate.QVERT %all equal within accuracy
        Qv = mean(Qv);
    else
        fprintf('Warning: inconsitency in vertical component of Scan %s  during call to Qvert.\n', scan.FILE);
    end
end


if QSign == -1 
    Qv = -Qv;
    %(as it has before been calculated for Q=kf-ki)
end