function erg = nxsgetvar(d,varname)

erg = [];

ind = find(strcmpi(d.data_scan.scanned_variables.variables_names.label,varname));
if ~isempty(ind)
    erg = d.data_scan.scanned_variables.data(:,ind(1));
end
