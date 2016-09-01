lat.a = 4.46;               % Give lattice constants, angles, and orientation
lat.b = 4.46;
lat.c = 4.46;
lat.alpha = pi/180 * 90;    
lat.beta  = pi/180 * 90;
lat.gamma = pi/180 * 90;
UB = UBmatrix(lat,[1,1,0],[0,0,1]);

kf = 4.1;

%%

scandef = 'bs qh -1 -1 0 -20 dqh 0 0 0 .5 np 50 ';  % Give the scan definition

hx = [10,0,0];   % Give the field components (as in "dr hx .. .. ..")


%%

limit = 1;   % Limit for low current (1 corresponds to 1A for k=2.662)
Imax = 10;   % Maximum acceptable current in A


%% Do the calculation and plot...

figure
hold on

a3b = -60:60;  % a3b-range to consider
fieldmag = sqrt(sum(hx.^2));
fieldang = -acos(hx(1)/fieldmag);

clear lowi1ok highi1ok lowi3ok highi3ok kiok kfok

hkle = scanQs(scandef);  
[hbar, mass_n, meVJ] = getoption('hbar', 'mass_n', 'meVJ');
ki = sqrt(2*mass_n*hkle(:,4) / hbar^2 / meVJ + kf^2 * 1E20) * 1E-10;
[qx,qy,qz] = calcQS ( hkle(:,1), hkle(:,2), hkle(:,3), UB);
[a3,gu,gl,a3p,a4] = spectroangles(qx,qy,qz,ki,kf,16,1,[],'a3p=0');
q2 = ki.^2 + kf.^2 - 2*ki*kf.*cos(pi/180*a4');

sinbeta = kf ./ sqrt(q2') .* sin(pi/180*a4);
cosbeta = (-ki' + kf .* cos(pi/180*a4)) ./ sqrt(q2') ;
beta = -acos(cosbeta') * sign(sinbeta) + fieldang;


fprintf('Scan from HKLE = %d,%d,%d,%d  to  HKLE = %d,%d,%d,%d.\n', hkle(1,1), hkle(1,2), hkle(1,3), hkle(1,4), hkle(end,1), hkle(end,2), hkle(end,3), hkle(end,4));

%fprintf('  np    H      K      L      En     a4    \n');

for i = 1:size(hkle,1)
    
    Imin_f = limit .* kf / 2.662;
    Imin_i = limit .* ki(i) / 2.662;
    
    % Test if a3b collides with ki
    kicoll = ~ (abs(a3b) < 43);
    % Test if a3b collides with kf
    kfcoll = ~ ((a4(i) - a3b < 100) & (a4(i) - a3b > 20));
    % Calculate i1 and test if ok
    i1 = -6.68 * cos(a3b*pi/180 + beta(i)) * fieldmag/10; 
    i1toolow  = abs(i1) < Imin_i;
    i1toohigh = abs(i1) > Imax;
    % Calculate i3 and test if ok
    i3 =  6.68 * cos(a3b*pi/180 + beta(i) + pi/3)  * fieldmag/10;
    i3toolow  = abs(i3) < Imin_f;
    i3toohigh = abs(i3) > Imax;

    a3bok = ~ (kicoll | kfcoll | i1toolow | i1toohigh | i3toolow | i3toohigh);
    
    a3bmatrix(i,:) = a3bok(:)';
    
    % Construct boundary lines
    kiok(i,:) = [43, -43];
    kfok(i,:) = [a4(i)-20,  a4(i)-100];
    lowi1ok(i,:)  = 180/pi * ([-acos(Imin_i/6.68*10/fieldmag), acos(Imin_i/6.68*10/fieldmag), acos(Imin_i/6.68*10/fieldmag)-pi, pi-acos(Imin_i/6.68*10/fieldmag)] - beta(i));
    highi1ok(i,:) = 180/pi * ([-acos(Imax  /6.68*10/fieldmag), acos(Imax  /6.68*10/fieldmag), acos(Imax  /6.68*10/fieldmag)-pi, pi-acos(Imax  /6.68*10/fieldmag)] - beta(i));
    lowi3ok(i,:)  = 180/pi * ([-acos(Imin_f/6.68*10/fieldmag), acos(Imin_f/6.68*10/fieldmag), acos(Imin_f/6.68*10/fieldmag)-pi, pi-acos(Imin_f/6.68*10/fieldmag)] - beta(i)) - 60 ;
    highi3ok(i,:) = 180/pi * ([-acos(Imax  /6.68*10/fieldmag), acos(Imax  /6.68*10/fieldmag), acos(Imax  /6.68*10/fieldmag)-pi, pi-acos(Imax  /6.68*10/fieldmag)] - beta(i)) - 60 ;
    
%    fprintf('%4d %6.2f %6.2f %6.2f %6.2f %6.2f\n', i, hkle(i,1), hkle(i,2), hkle(i,3), hkle(i,4), a4(i) );
    if isnan(a4(i)), fprintf('Step %d impossible\n',i); end
end

%%

lowi1ok(imag(lowi1ok)~=0) = NaN;
lowi3ok(imag(lowi3ok)~=0) = NaN;
highi1ok(imag(highi1ok)~=0) = NaN;
highi3ok(imag(highi3ok)~=0) = NaN;

%pcolor(1:size(hkle,1),a3b,a3bmatrix'*20);
p1 = plot(kiok,'-b');
p2 = plot(kfok,'-k');
p3 = plot(lowi1ok,'-r');
p4 = plot(highi1ok,'-.r');
p5 = plot(lowi3ok,'-g');
p6 = plot(highi3ok,'-.g');

p7 = 0;
for i=1:size(hkle,1)
   if any(a3bmatrix(i,:)), p7 = plot(ones(1,sum(a3bmatrix(i,:)))*i, a3b(a3bmatrix(i,:)), '.r'); end
end

title(scandef);
legend([p1(1),p2(1),p3(1),p4(1),p5(1),p6(1),p7], 'Collision with ki',  'Collision with kf', 'I1 too low', 'I1 too high', 'I3 too low', 'I3 too high',  'Possible a3b');
xlabel('Scan point number');
ylabel('a3b');
ylim([-50,50]);

