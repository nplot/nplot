% Lattice information
lat.a = 6.141;
lat.b = 10.736;
lat.c = 5.986;
lat.alpha = 82.27 * pi/180;
lat.beta  =107.43 * pi/180;
lat.gamma =102.67 * pi/180;

% kf in Ang^-1
kf = 3; 

% Give here the two orienting vector ([ax,ay,az], [bx,by,bz]) to calculate the UB matrix
UB=UBmatrix(lat,[1,0,0],[0,1,0]);

%% Calculate the spectrometer angles for given H,K,L and E

% define H, K, L, E
HKLE=[2,0,0.1,0];   


ki = sqrt(.4826*HKLE(4)+kf^2); 
[qx,qy,qz] = calcQS(HKLE(1),HKLE(2),HKLE(3),UB); 
[a3,gu,gl,a3p,a4,gfc] = spectroangles(qx,qy,qz,ki,kf,31,+1,[],'A3P=0');
a4 = a4 + 37.5; %(because just calculated for channel 31)

fprintf('a3 = %f,  gu = %f,  gl = %f,  a4 = %f,  gfc = %f\n',a3,gu,gl,a4,gfc);

%%

% Backward calculation: give a3,gu,gl,a4

[qlx,qly,qlz] = QLab(a4,0,16,ki,kf);
[qx,qy,qz] = QSampleA3(a3,gu,gl,0,qlx,qly,qlz);

[H,K,L] = calcHKL(qx,qy,qz,UB);
disp([H,K,L]);



