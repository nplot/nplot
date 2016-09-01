function [val,asgn] = mergelist(datalist,field)

% single list for all entries of datalist (if cell array)
% field: name of field (e.g. 'valuelist' etc.)
% asgn : index of slice for each row
% (simpler than cmbavg)

% P. Steffens, 08/2014

val = []; asgn =[];
ndata = length(datalist);
if isempty(datalist), return; end
if ndata==1 && ~iscell(datalist), m=datalist; clear datalist; datalist{1}=m; clear m; end % if not, make cell array

for nd = 1:ndata
    if ~isfield(datalist{nd},field), val = []; return; end
    val = [val; datalist{nd}.(field)]; %#ok<*AGROW>
    asgn = [asgn; nd*ones(size(datalist{nd}.(field),1),1)]; 
end
