function calcspurion(scanfile,command)


% lat.a = 9.01;               % Give lattice constants, angles, and orientation
% lat.b = 10.3;
% lat.c = 9.08;
% lat.alpha = pi/180 * 90;    
% lat.beta  = pi/180 * 90;
% lat.gamma = pi/180 * 90;
% avec = [0,1,0];             % first orienting vector
% bvec = [0,0,1];             % second orienting vector
% 
% kf = 2.662;
% 
% %%
% 
% scandef = 'bs qh 0 3 0.2 4 dqh 0 0 0  .5 np 30 ';  % Give the scan definition

% hkle = scanQs(scandef);  
% [hbar, mass_n, meVJ] = getoption('hbar', 'mass_n', 'meVJ');
% ki = sqrt(2*mass_n*hkle(:,4) / hbar^2 / meVJ + kf^2 * 1E20) * 1E-10;

%% Calculate the scan and plot...

scan = tasread(scanfile);

% ** check if Q-scan ??

ki = getvar(scan,'KI');
kf = getvar(scan,'KF');
lat = getlattice(scan);
avec = [scan.PARAM.AX, scan.PARAM.AY, scan.PARAM.AZ];
bvec = [scan.PARAM.BX, scan.PARAM.BY, scan.PARAM.BZ];
UB = UBmatrix(lat,avec,bvec);
pointnum = scan.DATA.PNT;


[qx,qy,qz] = calcQS ( getvar(scan,'QH'), getvar(scan,'QK'), getvar(scan,'QL'), UB);
[a3,gu,gl,a3p,a4] = spectroangles(qx,qy,qz,ki,kf,16,1,[],'a3p=0');

q2 = ki.^2 + kf.^2 - 2*ki.*kf.*cosd(a4');

figure
annotation('textbox',[0,.95,1,.05],'string',scan.COMND,'horizontalalignment','center','linestyle','none','fontsize',15);


subplot(3,3,3)
plot(pointnum,sqrt(q2),'-o'); ylabel('Qmod (A^{-1})'); yl=ylim; ylim([yl(1)-.1,yl(2)+.1]);
subplot(3,3,6)
plot(pointnum,a4,'-o'); ylabel('a4'); grid on
subplot(3,3,9)
plot(pointnum,ki,'-o'); ylabel('ki'); 
xlabel('Scan point number'); grid on


subplot(3,3,[1,8]);
hold on
xlabel('Qx (A^{-1})'); ylabel('Qy (A^{-1})');


zmax = max(qz); zmin = min(qz);
plot3(qx,qy,round(qz*1E7)*1E-7,'ob','Tag','Nominal scan points');

% Process kf' = KI (inc. on analyzer, elastic on sample)
QLx = -ki.* sind(a4');
QLy = ki .*(cosd(a4')-1);
QLz = zeros(size(QLx));
[qx,qy,qz] = QSampleA3(a3,gu,gl,a3p,QLx,QLy,QLz);

zmax = max([zmax, max(qz)]); zmin = min([zmin, min(qz)]);
plot3(qx,qy,round(qz*1E7)*1E-7,'-or','Tag','elastic scattering with KI, incoh on Ana');


% Process kf' = 2KI (second order incoming beam, inc. on analyzer, elastic on sample)
QLx = -2*ki.* sind(a4');
QLy = 2*ki .*(cosd(a4')-1);
QLz = zeros(size(QLx));
[qx,qy,qz] = QSampleA3(a3,gu,gl,a3p,QLx,QLy,QLz);

zmax = max([zmax, max(qz)]); zmin = min([zmin, min(qz)]);
plot3(qx,qy,round(qz*1E7)*1E-7,'-oc','Tag','second order ki, elastic scattering, incoh on Ana');


% Process ki' = KF (inc. on mono, elastic on sample)
QLx = -kf.* sind(a4');
QLy = kf .*(cosd(a4')-1);
QLz = zeros(size(QLx));
[qx,qy,qz] = QSampleA3(a3,gu,gl,a3p,QLx,QLy,QLz);

zmax = max([zmax, max(qz)]); zmin = min([zmin, min(qz)]);
plot3(qx,qy,round(qz*1E7)*1E-7,'-og','tag','incoh. on Mono, elastic scattering with KF');

%% Draw orienting vectors
end1 = avec * UB';
end2 = bvec * UB';
quiver([0,0], [0,0], [end1(1),end2(1)], [end1(2),end2(2)], 0, 'k');
text(end2(1)+.1*end1(1),end2(2), num2str(bvec,'[%g,%g,%g]'),'Fontname','Verdana','clipping','on'); 
if sign(end1(1))==-1; al='right'; else al='left'; end
text(1.1*end1(1),end1(2),  num2str(avec,'[%g,%g,%g]'), 'Fontname','Verdana','horizontalalignment',al,'clipping','on');

%% Draw reciprocal lattice

axis equal
xl=xlim; yl=ylim; zl=zlim;
corner = [  xl(1), xl(1), xl(1), xl(1), xl(2), xl(2), xl(2), xl(2); ...
            yl(1), yl(1), yl(2), yl(2), yl(1), yl(1), yl(2), yl(2); ...
            zmin , zmax , zmin , zmax , zmin , zmax , zmin , zmax ];
hklc = UB \ corner;
[hh,kk,ll] = ndgrid(ceil(min(hklc(1,:))-1e-6):floor(max(hklc(1,:))+1e-6), ceil(min(hklc(2,:))-1e-6):floor(max(hklc(2,:))+1e-6), ceil(min(hklc(3,:))-1e-6):floor(max(hklc(3,:))+1e-6));
qreclat = UB*[hh(:)';kk(:)';ll(:)'];
inrange = qreclat(1,:)>=xl(1) & qreclat(1,:)<=xl(2) & qreclat(2,:)>=yl(1) & qreclat(2,:)<=yl(2) & qreclat(3,:)>=zmin-1e-6 & qreclat(3,:)<=zmax+1e-6;
qreclat = qreclat(:,inrange); hh = hh(inrange); kk=kk(inrange); ll=ll(inrange);
p=plot3(qreclat(1,:),qreclat(2,:),round(qreclat(3,:)*1e7)*1e-7,'.k','markersize',22,'color',[.5,.5,.5],'tag','Reciprocal lattice'); 
uistack(p,'bottom');
if size(qreclat,2)<=100 % write hkl labels
    for r=1:size(qreclat,2)
        text(qreclat(1,r),qreclat(2,r),round(qreclat(3,r)*1e7)*1e-7, num2str([hh(r),kk(r),ll(r)],'[%g,%g,%g]'), 'Fontname','Verdana','fontsize',8,'horizontalalignment','center','verticalalignment','bottom','color',[.7,.7,.7],'clipping','on');
    end
end
xlim(xl); ylim(yl); 
box on

%% Powder lines Aluminium

AlPowQ =  2*pi / 4.0497 * [sqrt(3), sqrt(4), sqrt(8), sqrt(11), sqrt(12), 4, sqrt(19), sqrt(20), sqrt(24), sqrt(27)]; 
AlPowLab = {'111','200','220','113','222','400','331','420','422','115'};

ang=0:.5:360;
for p=1:numel(AlPowQ)
    plot(AlPowQ(p)*cosd(ang), AlPowQ(p)*sind(ang), ':k', 'Tag', ['Al ', AlPowLab{p},' powder line']);
end



