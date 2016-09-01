function lat=getlattice(scan)

% Extract lattice information from MAD scan structure
% Use reciprocal to obtain reciprocal lattice parameters
%
% P. Steffens, 01/2008


if ischar(scan), scan = tasread(scan); end
scan = scan(1);

if all(fieldcheck(scan.PARAM, {'AS', 'BS', 'CS', 'AA', 'BB', 'CC'}))
    lat.a = scan.PARAM.AS;
    lat.b = scan.PARAM.BS;
    lat.c = scan.PARAM.CS;
    lat.alpha = pi/180 * scan.PARAM.AA;
    lat.beta  = pi/180 * scan.PARAM.BB;
    lat.gamma = pi/180 * scan.PARAM.CC;
else
    lat=[];
    disp('Structure does not contain all necessary fields.');
end
