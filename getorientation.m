function [avec,bvec]=getorientation(scan)

% Extract crystal orientation information from MAD scan structure
% Can be filename. If multiple, take first.
%
% P. Steffens, 03/2008


if ischar(scan), scan = tasread(scan); end
scan = scan(1);

if all(fieldcheck(scan.PARAM, {'AX', 'AY', 'AZ', 'BX', 'BY', 'BZ'}))
    avec = [scan.PARAM.AX, scan.PARAM.AY, scan.PARAM.AZ];
    bvec = [scan.PARAM.BX, scan.PARAM.BY, scan.PARAM.BZ];
else
    avec=[];
    bvec=[];
    disp('Structure does not contain all necessary fields.');
end
