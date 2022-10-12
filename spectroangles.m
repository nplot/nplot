function [a3,gu,gl,a3p,a4,gfc]=spectroangles(vx,vy,vz,ki,kf,ch,SS,mintilt,varargin)

% Calculate spectrometer angles to bring Q-vector in detector No. ch
%
% Q=(vx,vy,vz) in crystal cartesian coordinates UB*(hkl)
% vx..kf as usual; give "ch" only once (single value, no vector)!
% varargin (optional) may contain:
%   opt:      desired solution (give condition, e.g. 'A3P=0', 'A3=25.4' etc.)
%   zerovals: angles zeros
%
% P. Steffens, 11/2008 - 1/2015

QSign = getoption('QSign');
if isempty(SS), SS=getoption('stdSS'); end
if isempty(mintilt), mintilt = false;  end

if QSign==-1
    vx=-vx; vy=-vy; vz=-vz;
    %if Q=ki-kf, invert coord. because following calculation is for Q=kf-ki
end

[ok,vx,vy,vz,ki,kf] = makesamesize(vx,vy,vz,ki,kf);
if ~ok, return; end

nu  = pi/180*(15+(2.5*(ch-1)));
b   = (ki.^2+kf.^2-(vx.^2+vy.^2+vz.^2))./(2*ki.*kf);
c   = vz ./ ki;

ind = find(((1-b.^2-c.^2)<0)+(b<=-sin(nu)*sqrt(1-c.^2)));
%These are not accessible!
a4(ind) = NaN; gfc(ind)= NaN;
a3(ind) = NaN; a3p(ind)= NaN;
gu(ind) = NaN; gl(ind) = NaN;

%Continue with the rest
ind = find(((1-b.^2-c.^2)>=0).*(b>-sin(nu).*sqrt(1-c.^2)));
vx  = vx(ind);   vy = vy(ind);   vz = vz(ind);
ki  = ki(ind);   kf = kf(ind); 
b   = b(ind);     c = c(ind);

nu    = pi/180*(15+(2.5*(ch-1)));
cosnu = cos(nu);
sinnu = sin(nu);


if SS == 1 
    sindelta = - b * sinnu + cosnu * sqrt(1 - c.^2 - b.^2);
    delta = asin(sindelta);
else
    sindelta = - b * sinnu - cosnu * sqrt(1 - c.^2 - b.^2);
    % delta may be larger or smaller than -90. Regard sign of cos(delta):
    
    % correction:
    posind = b + sindelta*sinnu > 0;
    delta = 0*sindelta; % (initialize with same size)
    delta(posind) = asin(sindelta(posind));
    delta(~posind)= -pi - asin(sindelta(~posind));
    
%     if b + sindelta*sinnu > 0
%         delta = asin(sindelta);
%     else
%         delta = -pi - asin(sindelta);
%     end
end

cosdelta  = cos(delta);
singamma  = c ./ cosdelta;
cosgamma  = sqrt(1 - singamma.^2);

a4(ind)   = 180/pi*delta+52.5;


gfc(ind)  = 180/pi*asin(singamma);

% if abs(b-cosdelta*cosgamma*cosnu+sindelta*sinnu)>1E-10 || abs(cosdelta*singamma-c) > 1E-10, disp('Fehler in Winkelberechnung!'); end %Test

j1mod = sqrt(-2*kf.*cosdelta.*cosgamma*cosnu.*ki + kf.^2+2*kf.*sindelta*sinnu.*ki + ki.^2-ki.^2.*cosdelta.^2+ki.^2.*cosgamma.^2.*cosdelta.^2);

for m=1:numel(ind) % for all possible, do the calculation
    
    % Construct two right-handed frames to obtain rotation matrix
    % (How, depends on if Flatcone or minimum-tilt solution is desired
    
    if mintilt
        i0 = [-vy(m);  vx(m); vz(m)] / sqrt(vx(m)^2+vy(m)^2+vz(m)^2);
        j0 = [-vx(m); -vy(m);     0] / sqrt(vx(m)^2+vy(m)^2);
        k0 = cross(i0,j0);   k0 = k0 / sqrt(k0'*k0);
        
        i1 = [-kf(m) * cosdelta(m) ; -kf(m) *sindelta(m) - ki(m); 0];   i1 = i1 / sqrt(i1'*i1);
        j1 = [-i1(2); i1(1); 0];
        k1 = cross(i1,j1);   k1 = k1 / sqrt(k1'*k1);        
    
    else
        
        k0 = [0; 0; 1];
        j0 = [-vx(m); -vy(m); 0] / sqrt(vx(m)^2+vy(m)^2);
        i0 = [-vy(m);  vx(m); 0] / sqrt(vx(m)^2+vy(m)^2);

        k1 = [sindelta(m)*singamma(m); -cosdelta(m)*singamma(m); cosgamma(m)];
        j1 = 1/j1mod(m) * [-cosdelta(m)*kf(m)*cosnu+cosgamma(m)*kf(m)*sindelta(m)*sinnu+cosgamma(m)*ki(m); ...
                           -kf(m)*(cosgamma(m)*cosdelta(m)*sinnu+sindelta(m)*cosnu); ...
                           -singamma(m)*(kf(m)*sinnu+sindelta(m)*ki(m))];
        i1 = 1/j1mod(m) * [-kf(m)*cosgamma(m)*sindelta(m)*cosnu-cosdelta(m)*kf(m)*sinnu-cosdelta(m)*sindelta(m)*ki(m)+cosdelta(m)*sindelta(m)*ki(m)*cosgamma(m)^2; ... 
                           -sindelta(m)*kf(m)*sinnu-ki(m)+ki(m)*cosdelta(m)^2-ki(m)*cosgamma(m)^2*cosdelta(m)^2+cosgamma(m)*cosdelta(m)*kf(m)*cosnu; ...
                           -singamma(m)*(cosdelta(m)*cosgamma(m)*ki(m)-kf(m)*cosnu)];
    end
       
    
    M = [i1, j1, k1] * ([i0, j0, k0]^(-1));  % This is the desired rotation
    
    % Now, let the subroutine do the rest
    if m==1
        [a3(ind(m)),gu(ind(m)),gl(ind(m)),a3p(ind(m)),opt,zerovals,config] = anglesfrommatrix(M, varargin{:});
    else
        [a3(ind(m)),gu(ind(m)),gl(ind(m)),a3p(ind(m))] = anglesfrommatrix(M,opt,zerovals,config);
    end
   
end
