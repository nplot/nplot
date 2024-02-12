function nxsprintscan(nxs,outdef)

% output scan data in text format
%
% by default, use "out" list from Nomad (saved in Nexus File)
% if "outdef" given, use this instead, e.g.: nxsprintscan(d,'A3 A4 TRT')

% PS 06/23 - 02/24

% scanned Vars
varlist = fieldnames(scansteps(nxs.(nxs.instrument_name).command_line.actual_command));

% counts
varlist = [varlist(:)', {'Monitor1'}, {'Monitor2'}, {'Time'}, {'Detector'}];

% special vars
if any(strcmpi(nxs.data_scan.scanned_variables.variables_names.label,'F1'))
    varlist = [varlist(:)', {'F1'}];
end
if any(strcmpi(nxs.data_scan.scanned_variables.variables_names.label,'F2'))
    varlist = [varlist(:)', {'F2'}];
end

% "out"
if nargin>1 && ~isempty(outdef)
    outlist = strsplit(outdef,{',',' '});
else
    outlist = strsplit(nxs.(nxs.instrument_name).command_line.output_list,',');
end
varlist = [varlist(:)', outlist(:)'];

% name conventions
varlist = replacename(varlist,'HX','HelmholtzHX');
varlist = replacename(varlist,'HY','HelmholtzHY');
varlist = replacename(varlist,'HZ','HelmholtzHZ');

% Extract values
mat = [];
for vi = 1:length(varlist)
    val = nxsgetvar(nxs,varlist{vi});
    if isempty(val)
        error(['Cannot extract column ' varlist{vi}]);
    end
    mat = [mat, val]; %#ok<AGROW>
end
mat = [(1:size(mat,1))', mat];

% add "pnt"
varlist = [{'PNT'}, varlist(:)'];
% replace names according to TAS convention
varlist = replacename(varlist,'Monitor1','M1');
varlist = replacename(varlist,'Monitor2','M2');
varlist = replacename(varlist,'Time','TIME');
varlist = replacename(varlist,'Detector','CNTS');
varlist = replacename(varlist,'HelmholtzHX','HX');
varlist = replacename(varlist,'HelmholtzHY','HY');
varlist = replacename(varlist,'HelmholtzHZ','HZ');


% Print column headers
for vi =1:length(varlist)
    fprintf('%12s', varlist{vi});
end
fprintf('\n');

% Print data
for zi = 1:size(mat,1)
    for vi =1:size(mat,2)
        if any(strcmpi(varlist{vi},{'PNT','M1','M2','CNTS','F1','F2'}))
            fprintf('  %10d', mat(zi,vi));
        elseif strcmpi(varlist{vi},'TIME') % various precisions
            fprintf('    %8.2f', mat(zi,vi));
        elseif nxs.data_scan.scanned_variables.variables_names.axis(strcmpi(nxs.data_scan.scanned_variables.variables_names.label,varlist{vi}))
            fprintf('    %8.3f', mat(zi,vi));
        else
            fprintf('    %8.4f', mat(zi,vi));
        end
    end
    fprintf('\n');
end



function varlist = replacename(varlist,old,new)
for ii=find(strcmpi(varlist,old))
    varlist{ii} = new;
end






