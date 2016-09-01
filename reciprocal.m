function rezlat=reciprocal(lat)

% Do the transform to reciprocal lattice
% (angles in radian)

% P. Steffens, 01/2008

V=lat.a*lat.b*lat.c*sqrt(1-cos(lat.alpha)^2-cos(lat.beta)^2-cos(lat.gamma)^2+2*cos(lat.alpha)*cos(lat.beta)*cos(lat.gamma));

rezlat.a=2*pi/V*lat.b*lat.c*sin(lat.alpha);
rezlat.b=2*pi/V*lat.a*lat.c*sin(lat.beta);
rezlat.c=2*pi/V*lat.a*lat.b*sin(lat.gamma);

rezlat.alpha = acos((cos(lat.beta) *cos(lat.gamma)-cos(lat.alpha))/(sin(lat.beta) *sin(lat.gamma)));
rezlat.beta  = acos((cos(lat.alpha)*cos(lat.gamma)-cos(lat.beta)) /(sin(lat.alpha)*sin(lat.gamma)));
rezlat.gamma = acos((cos(lat.alpha)*cos(lat.beta) -cos(lat.gamma))/(sin(lat.alpha)*sin(lat.beta )));

rezlat.V0 = sqrt(1-cos(rezlat.alpha)^2-cos(rezlat.beta)^2-cos(rezlat.gamma)^2+2*cos(rezlat.alpha)*cos(rezlat.beta)*cos(rezlat.gamma));
rezlat.V = rezlat.a*rezlat.b*rezlat.c*rezlat.V0;
