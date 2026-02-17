% Script to simulate a simple Flatcone Scan


sample.lattice.a = 5.1;
sample.lattice.b = 3.8;
sample.lattice.c = 6.2;
sample.lattice.alpha = 90 * pi/180;
sample.lattice.beta  = 90 * pi/180;
sample.lattice.gamma = 90 * pi/180;
sample.ax = [1 0 0]; % first vector defining scattering plane
sample.bx = [0 1 0]; % second
    
a3 = 10:1:50;   % A3- range
a4 = -50;       % A4 (fixed)
EN = 1;         % Energy transfer

clear d
[a3,a4] = ndgrid(a3, a4 + linspace(-37.5,37.5,31));
d.coordlist = [a4(:), a3(:)];
d.valuelist = 100 * rand(size(d.coordlist(:))); d.valuelist(:,2) = sqrt(d.valuelist);
d.type = 'CONST-ENERGY, CONST-QVERT';
d.coordtype = 'ANGLES';
d.raw = true;
d.sampleinfo = sample;
d.KF = 1.4;
d.KI = sqrt((2.0721*d.KF^2 + EN)/2.0721);
d.QVERT = 0;
d.constants = {'KF'  'KI'  'QVERT'};
d.expname = 'Flatcone Scan';
d.dataname = 'Simulation';
d.monitorlist = ones(size(d.valuelist));

fcplot(d);


