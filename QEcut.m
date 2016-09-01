function [H,K,L,en,data] = QEcut(figurelist,scandef)

% Extract data from open figure windows (give numbers in "figurelist")
% that lie along a scan line defined by "scandef" (qvert and energy of
% scandef are ignored)


HKLE = scanQs(scandef);
H = HKLE(:,1); K = HKLE(:,2); L = HKLE(:,3);
d1 = getfiguredata(figurelist(1));
UB = UBmatrix( d1.sampleinfo.lattice, d1.sampleinfo.ax, d1.sampleinfo.bx);
[qx,qy] = calcQS(H,K,L,UB);
[hbar, mass_n, meVJ] = getoption('hbar', 'mass_n', 'meVJ');

%% Get data from interpolation in each map

for nfig = 1:numel(figurelist)
    dat = getfiguredata(figurelist(nfig));
    en(nfig) = (dat.KI^2-dat.KF^2)*1E20*hbar^2/2/mass_n*meVJ; % E-transfer
    datxy = coordtransform(dat,'qxy');
    data(nfig,:) = linearinterpolation(datxy,[qx,qy]);
end
[en,order] = sort(en);
data = data(order,:);

%% Plot the data

xvar = H;  % Which variable for x-axis?

xvert = (xvar([1,1:end])+xvar([1:end,end]))/2;
envert= (en([1,1:end])+en([1:end,end]))/2;
figure
pcolor(xvert,envert,log([data,data(1:end,end);data(end,1:end),data(end,end)]));
shading flat
ylabel('Energy (meV)');