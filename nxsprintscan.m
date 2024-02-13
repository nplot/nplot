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
    outlist = strsplit(strtrim(outdef),{',',' '});
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
        val = nxsgetvar(nxs,varlist{vi},'calc'); % try again, extended search
        if isempty(val)
            fprintf(2,['Cannot extract column ' varlist{vi} '\n']);
            return;
        else
            warning(['Column ' varlist{vi} ' not in scan data, but obtained indirectly.']);
        end
    end
    mat = [mat, val]; %#ok<AGROW>
end
% Add PNT (and PAL)
ndat = size(mat,1);
if hasfield(nxs.(nxs.instrument_name),'pal') && hasfield(nxs.(nxs.instrument_name).pal,'pal_contents'), npal = numel(strfind(upper(nxs.(nxs.instrument_name).pal.pal_contents),'CO ')); else npal=0; end
if npal >0 % polarized
    % add "pal"
    varlist = [{'PAL'}, varlist(:)'];  % ** ! is "guessed". Improve when PAL becomes available in nexus!!
    mm = repmat((1:npal)',1,ceil(ndat/npal));
    mat = [mm(1:ndat)', mat];
else 
    npal = 1;
end
% add "pnt"
varlist = [{'PNT'}, varlist(:)'];
mm = repmat(1:ndat,npal,1);
mat = [mm(1:ndat)', mat];

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
    if any(strcmpi(varlist{vi},{'PNT','PAL','F1','F2'}))
        fprintf('%6s', varlist{vi});
    else
        fprintf('%12s', varlist{vi});
    end
end
fprintf('\n');

% Print data
for zi = 1:size(mat,1)
    for vi =1:size(mat,2)
        if any(strcmpi(varlist{vi},{'PNT','PAL','F1','F2'}))
            fprintf(' %5d', mat(zi,vi));
        elseif any(strcmpi(varlist{vi},{'M1','M2','CNTS'}))
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






