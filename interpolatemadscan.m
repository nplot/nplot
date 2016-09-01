function [z,dz] = interpolatemadscan(datalist, scancommand)

% [z,dz] = interpolatemadscan(datalist, scancommand)
%
% Obtain an interpolated scan by providing a scan-command in the MAD syntax.
% So far, works if datalist is a constant energy cut.
%
% Example: [z,dz]=interpolatemadscan(getfiguredata, 'sc qh 1 0 0 0 dqh 0 .01 0 0 np 21');

% P. Steffens, 07/2010


Qpos = scanQs(scancommand);

datalist = coordtransform(datalist,'qxy');

UB = UBmatrix( datalist.sampleinfo.lattice, datalist.sampleinfo.ax, datalist.sampleinfo.bx );
for i=1:size(Qpos,1)
    Qxypos(i,1:3) = (UB * Qpos(i,1:3)')';
end

[mass_n, meVJ, hbar, maxdeviate] = getoption('mass_n','meVJ','hbar','maxdeviate');

if isfield(datalist,'QVERT')    % check for vertical Q
    if max(abs(datalist.QVERT - Qxypos(:,3))) > maxdeviate.QVERT
        fprintf('Warning: The vertical component of Q does not match. Please check the result!!\n(You may increase maxdeviate.QVERT in options.m to supress this warning)\n');
    end
end

if isfield(datalist,'KI') && isfield(datalist,'KF')     %check energy
    en = hbar^2/2/mass_n * 1E20 * meVJ * (datalist.KI^2-datalist.KF^2);
    if max(abs(en - Qpos(:,4))) > hbar^2/2/mass_n * (maxdeviate.KI^2 + maxdeviate.KF^2)  * 1E20 * meVJ 
        fprintf('Warning: The deviation in energy seems too high. Please check the result!!\n(You may increase maxdeviate.KI/KF in options.m to supress this warning)\n');
    end
end

[z, dz] = linearinterpolation(datalist, Qxypos(:,1:2));

