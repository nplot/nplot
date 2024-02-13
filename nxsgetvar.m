function erg = nxsgetvar(d,varname,varargin)

% obtain scan data from nexus file
% 
% nxsgetvar(.., .., 'calc') tries to deduce, if not present in scan data

% PS 02/24

erg = [];

ind = find(strcmpi(d.data_scan.scanned_variables.variables_names.label,varname));
if ~isempty(ind)
    erg = d.data_scan.scanned_variables.data(:,ind(1));
elseif any(strcmpi(varargin, 'calc')) % try to obtain by searching and calculation
    erg = getvar(nxsTAS(d),varname,'includeparam','readheader');    
end
