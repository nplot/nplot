function erg = nxsgetvar(d,varname,varargin)

% erg = nxsgetvar(d,varname,varargin)
% 
% Obtain scan data from nexus file
% 
% nxsgetvar(.., 'A1.target_value') looks in subfields
% nxsgetvar(.., .., 'calc') tries to deduce, if not present in scan data

% PS 02/26

erg = [];
ind = [];

if isfield(d,'data_scan') && isfield(d.data_scan,'scanned_variables') && isfield(d.data_scan.scanned_variables,'variables_names') && isfield(d.data_scan.scanned_variables.variables_names,'label')
    ind = find(strcmpi(d.data_scan.scanned_variables.variables_names.label,varname));   % varname is in scan ?
end

if ~isempty(ind)
    erg = d.data_scan.scanned_variables.data(:,ind(1));

else    % otherwise, look in param list
    s = strsplit(varname,'.');
    varname = s{1};
    if length(s)>1, subfield = s{2}; else, subfield = 'value'; end
    if upper(varname(1))=="Z", varname = varname(2:end); subfield = 'offset_value';
    else, switch upper(varname)
        case 'DM', varname='Monochromator'; subfield='d_spacing'; 
        case 'DA', varname='Analyser'; subfield='d_spacing'; 
        case {'EI','KI'}, varname='Monochromator'; subfield=lower(varname); 
        case {'EF','KF'}, varname='Analyser'; subfield=lower(varname); 
    end, end
    paramlist = fieldnames(d.(d.instrument_name));
    if any(strcmpi(paramlist,varname))
        paramname = paramlist{(strcmpi(paramlist,varname))};% this is now correct upper/lower case
    
        if isfield(d.(d.instrument_name),paramname) && isfield(d.(d.instrument_name).(paramname),subfield)
            erg = d.(d.instrument_name).(paramname).(subfield);
        elseif isfield(d.(d.instrument_name),paramname)
            erg = d.(d.instrument_name).(paramname) ;
        end
    elseif any(strcmpi(varargin, 'calc')) % try to obtain by searching and calculation
        erg = getvar(nxsTAS(d),varname,'includeparam','readheader');  
    end

% elseif isfield(d.(d.instrument_name),upper(varname))
%     erg = d.(d.instrument_name).(upper(varname));
% elseif any(strcmpi(varname,{'ki','ei'})) && isfield(d.(d.instrument_name),'Monochromator')
%     erg = d.(d.instrument_name).Monochromator.(lower(varname));
% elseif any(strcmpi(varname,{'kf','ef'})) && isfield(d.(d.instrument_name),'Analyser')
%     erg = d.(d.instrument_name).Analyser.(lower(varname));
  
end
