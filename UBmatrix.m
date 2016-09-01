function [UB,U,B]=UBmatrix(lat,ref1,ref2)

%Calculate the UB matrix
%lat:       crystal lattice information (a,b,c,alpha,beta,gamma)
%ref1,ref2: h,k,l of two reciprocal lattice vectors in the scattering plane
%
%note that ref1 // ki for A3=0 
%
% P. Steffens, 07/2008

rezlat=reciprocal(lat);

ref1=ref1(:);
ref2=ref2(:); %ensure column vectors

%reciprocal lattice vectors in the crystal cartesian system
astar = rezlat.a * [1;0;0];
bstar = rezlat.b * [cos(rezlat.gamma); sin(rezlat.gamma); 0];
cstar = rezlat.c * [cos(rezlat.beta);  (cos(rezlat.alpha)-cos(rezlat.beta)*cos(rezlat.gamma))/sin(rezlat.gamma); rezlat.V0/sin(rezlat.gamma)];

%now follow the concept of Busing & Levy, 1967

%unit vector of first orienting reflection in crystal cartesian coordinates
t1 = ref1(1)*astar + ref1(2)*bstar + ref1(3)*cstar;
t1 = t1 / sqrt(t1'*t1);
%second orthogonal unit vector in the scattering plane
t2 = ref2(1)*astar + ref2(2)*bstar + ref2(3)*cstar;
t2 = t2 - t1 * (t1'*t2);
t2 = t2 / sqrt(t2'*t2);
%complete right-hand orthogonal system
t3 = cross(t1,t2);

Tc = [t1, t2, t3];

U = Tc'; %=Tc^(-1)

B = [ rezlat.a,  rezlat.b * cos(rezlat.gamma),  rezlat.c * cos(rezlat.beta); ...
      0,         rezlat.b * sin(rezlat.gamma), -rezlat.c * sin(rezlat.beta)*cos(lat.alpha); ...
      0,         0,                             2*pi/lat.c ];
  
UB = U*B;
