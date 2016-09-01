function time = totaltime(scans)

% calculate total counting time of scans 
% (** not yet for fcu data files)

[scans,nscans] = tasread(scans,'download');

time = 0;
for i=1:nscans
    if isstruct(scans), scan = scans(i); else scan = scanst{i}; end
    if isfield(scan,'DATA') && isfield(scan.DATA,'TIME')
        time = time + sum(scan.DATA.TIME);
    end
end

fprintf('Total counting time: %d h %d min %d sec \n', floor(time/3600), floor(mod(time,3600)/60), floor(mod(time,60)));
