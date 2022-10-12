function output=cmbavg(scanlist,opt,varargin)

% function output=cmbavg(scanlist,opt,varargin)
%
% Combine array of scans (linear list format) into a single one.
% This includes averaging unless opt="noAvg" is given
% opt="explicit": use grid given in varargin{1}
% opt="auto":     try to find optimum grid automatically (not yet implemented)
% opt="nogrid":   averaging without grid (ie. only exactly equal points)
% opt="standard" or empty: use standard grid parameters from options file
%
% varargin can contain a number of parameter pairs (override values from options.m)
% varargin containing 'ignore [CONST]' ignores deviations in constant(s)


% % P. Steffens 10/2016

if nargin<2 || isempty(opt),  opt = 'standard'; end 

output = [];
nscans = length(scanlist);
if isempty(scanlist), return; end
if nscans==1 && ~iscell(scanlist), m=scanlist; clear scanlist; scanlist{1}=m; clear m; end % if not, make cell array
if isempty(scanlist{1}), return; end

%Check if all lists are of the same type
type1 = upper(scanlist{1}.type);
coordtype1 = upper(scanlist{1}.coordtype);
raw1 = scanlist{1}.raw;
isnotpol1 = ~isfield(scanlist{1},'pallist') || (isfield(scanlist{1},'polarized') && scanlist{1}.polarized == false);
for i=2:nscans
    if ~strcmpi(scanlist{i}.type, type1)
        fprintf('Error: tried to combine data lists of different type! This is impossible. \n');
        return;
    end
    if ~strcmpi(scanlist{i}.coordtype, coordtype1)
        fprintf('Error: tried to combine data lists with different types of corrdinates! This is impossible. \n');
        return;
    end
    if scanlist{i}.raw ~= raw1
        fprintf('Error: tried to combine data lists of different value types (.raw = 0/1)! This is impossible. \n');
        return;
    end
    isnotpoli = ~isfield(scanlist{i},'pallist') || (isfield(scanlist{i},'polarized') && scanlist{i}.polarized == false);
    if isnotpol1 ~= isnotpoli
       fprintf('Error: tried to combine data lists with and without polarization! This is impossible. \n');
       return;
    end 
end
output.type = upper(type1);
output.coordtype = upper(coordtype1);
output.raw = raw1;
output.polarized = ~isnotpol1;

% expname and title
output.expname = scanlist{1}.expname;   output.dataname = scanlist{1}.dataname;
for i=1:nscans, enames{i}=scanlist{i}.expname; dnames{i}=scanlist{i}.dataname; end
enames = unique(enames,'stable'); dnames = unique(dnames,'stable');
if iscell(enames) && length(enames)>1, for i=2:length(enames), output.expname = [output.expname, ', ', enames{i}]; end; end
if iscell(dnames) && length(dnames)>1, for i=2:length(dnames), output.dataname= [output.dataname,', ', dnames{i}]; end; end

% properties...
propnames = {};
for i=1:nscans, if isfield(scanlist{i},'properties'), propnames = [propnames; fieldnames(scanlist{i}.properties)]; end, end
propnames = unique(propnames);
for p=1:length(propnames)
    pval={}; 
    for i=1:nscans
        if ~isfield(scanlist{i},'properties') || ~isfield(scanlist{i}.properties,propnames{p}), pval = {'not defined for all data'}; break; 
        else pval = [pval, scanlist{i}.properties.(propnames{p})]; end
    end
    if length(unique(pval))>1, output.properties.(propnames{p}) = 'not equal for all data';
    else pval = unique(pval);  output.properties.(propnames{p}) = pval{1}; end
end

        
        
if isfield(scanlist{1},'sampleinfo'), output.sampleinfo = scanlist{1}.sampleinfo; end

if isfield(scanlist{1},'paldeflist'), output.paldeflist = scanlist{1}.paldeflist; end
% ** Attention: no consitency check of the paldeflists is done! For
% polarized data, the points are just combined by their pal-values.


%Check consistency of constants
[eqset, eqvalues, deviate, notall] = checkconstants(scanlist);
%Create fields and assign values
for j=1:length(eqset),  output.(eqset{j}) = eqvalues{j};  end
output.constants = eqset;
%Warnings for missing values
for j=1:length(notall), fprintf('Warning: Value for "%s" does not exist in all data sets! Continue... (Please check!)\n', notall{j}); end
%Error in case of too large deviations
if ~isempty(deviate)
    errtxt=[]; 
    ignoreconsts = readinput('ignore',varargin); % names of constants to ignore (can be cell array)
    for j=1:length(deviate)
        if ~any(strcmpi(deviate{j},ignoreconsts)), errtxt = [errtxt, deviate{j}, ' ']; end
    end
    if ~isempty(errtxt)
        fprintf('Error while combining data lists: too large deviation in %s !\n (You may increase the maximum acceptance in the options file.) \n', errtxt);
        output=[];
        return;
    end
end


% Do the actual combination of the lists
coordlist = [];
valuelist = [];
vertexlist= [];
faces = [];
delaunaytri = [];
monitorlist = [];
sectionlist = [];
taglist={}; hastags = false;
pallist = [];
for i=1:nscans
     if (strcmpi(opt,'NOAVG')) 
        if isfield(scanlist{i},'faces'), 
            faces(size(faces,1)+(1:size(scanlist{i}.faces,1)),1:size(scanlist{i}.faces,2)) = scanlist{i}.faces + size(vertexlist,1); 
        end
        if isfield(scanlist{i},'delaunaytri'), delaunaytri = [delaunaytri; scanlist{i}.delaunaytri + size(coordlist,1)]; end %#ok<*AGROW>
        if isfield(scanlist{i},'vertexlist'), vertexlist = [vertexlist;scanlist{i}.vertexlist]; end    
        % Note that the combination of vertexlist, faces and delaunaytri does not necessarily make sense. 
        % It is intended for the case of several separate slices.
    end
    coordlist  = [coordlist; scanlist{i}.coordlist];
    valuelist  = [valuelist; scanlist{i}.valuelist];
    if output.raw, monitorlist = [monitorlist; scanlist{i}.monitorlist]; end
    if isfield(scanlist{i},'taglist'),
        taglist = [taglist, scanlist{i}.taglist];  % {taglist{:}, scanlist{i}.taglist{:}};
        hastags = true;
    else
        fill = cell(1,size(scanlist{i}.coordlist,1));
        taglist = [taglist, fill];  % fill with empty cells
    end
    if isfield(scanlist{i},'sectionlist')
        sectionlist = [sectionlist; scanlist{i}.sectionlist];
    end
    if output.polarized, pallist = [pallist; scanlist{i}.pallist]; end
end


% For psi-coordinate, bring to right range (there may be jumps)
if any(strcmpi(output.coordtype, {'ANGLES', 'QPLANE', 'ANGLESQZ'}))
    coordlist(:,2) = bringtorange(coordlist(:,2));
elseif any(strcmpi(output.coordtype, {'ANGLESENERGY', 'ENERGY3D'}))
    coordlist(:,3) = bringtorange(coordlist(:,3));
end

if (strcmpi(opt,'NOAVG'))  % (no averaging)
    output.coordlist = coordlist;
    output.valuelist = valuelist;
    if ~isempty(vertexlist), output.vertexlist = vertexlist; end
    if ~isempty(delaunaytri), output.delaunaytri = delaunaytri; end
    if ~isempty(faces), faces(faces==0)=NaN; output.faces = faces; end
    if ~isempty(sectionlist), output.sectionlist = sectionlist; end
    if output.raw, output.monitorlist = monitorlist; end
    if hastags, output.taglist = taglist; end
    if output.polarized, output.pallist = pallist; end
    return; 
end

% The combination is now done.
% Next, binning and averaging:

% Get some parameters either from varargin or options file 
bindist = readinput('bin',varargin);
if isempty(bindist) 
    stdbindist = getoption('stdbindist','check',varargin); 
    try bindist=stdbindist.(output.coordtype); 
    catch
        fprintf('Error: your options.m file does not contain the field ''stdbindist.%s''.\n',output.coordtype); output=[]; return;
    end
end
normval = readinput('monitor',varargin); if isempty(normval),normval = readinput('time',varargin); end
if isempty(normval), normval = getoption('normval','check',varargin); end
gridstep = readinput('grid',varargin);
%if isempty(gridstep) && ~strcmpi(opt,'EXPLICIT'), stdgrid = getoption('stdgrid','check',varargin); gridstep=stdgrid.(output.coordtype); end
% ** ??
if isempty(gridstep) 
    stdgrid = getoption('stdgrid','check',varargin); 
    try gridstep=stdgrid.(output.coordtype); 
    catch
        fprintf('Error: your options.m file does not contain the field ''stdgrid.%s''.\n',output.coordtype); output=[]; return;
    end
end


if strcmpi(output.coordtype,'GENERAL')
    if numel(bindist) < size(coordlist,2), bindist = bindist(1)* ones(1,size(coordlist,2)); end
    if numel(gridstep) < size(coordlist,2), gridstep = gridstep(1)* ones(1,size(coordlist,2)); end
end

% First, obtain a suitable grid for averaging the data

if output.polarized % Add pal-state as "coordinate"
    coordlist = [coordlist, pallist];
    bindist = [bindist, .1];
    gridstep = [gridstep, 1];
end

if strcmpi(opt,'EXPLICIT')
     grid = varargin{1};
     if output.polarized, pgrid = zeros(0,size(grid,2)); for p=unique(pallist)', pgrid = [pgrid; grid, p*ones(size(grid,1),1)]; end, grid=pgrid; end
     binlist = bintogrid(coordlist, grid, bindist);
elseif strcmpi(opt,'STANDARD')
%    [grid,assign] = makegrid(coordlist, median(coordlist, 1), gridstep);
    % ** Until better solution, take median as center and standard values as grid spacing
    [grid,assign] = makegrid(coordlist, [], gridstep);
    binlist = [(1:numel(assign))', assign];
elseif strcmpi(opt,'NOGRID')
    [grid, ~, assign] = unique(coordlist, 'rows');
    binlist = [(1:numel(assign))', assign];
elseif strcmpi(opt,'AUTO')
    fprintf('cmbavg: Option "auto" not yet implemented... Choose other option. \n');
    return;
end

% Finally, do the averaging. Points assigned to the same grid point are averaged
[ubl2, ~ , ubl2num] = unique(binlist(:,2));
output.coordlist = zeros(numel(ubl2),size(coordlist,2));
output.valuelist = zeros(numel(ubl2), 2);
if hastags, output.taglist = cell(numel(ubl2), 1); end

if output.raw ==1  % For counting data: just add, normalize and obtain error as sqrt  
    output.monitorlist = zeros(size(output.valuelist,1),size(monitorlist,2)); 
    output.coordlist = grid( ubl2, : );
    valuelist = valuelist(:,1) .* monitorlist(:,1) / normval;  % These are the real original counts!
    for i=1:numel(binlist(:,1))
        output.valuelist  (ubl2num(i),1) = output.valuelist(ubl2num(i))     + valuelist(i);
        output.monitorlist(ubl2num(i),:) = output.monitorlist(ubl2num(i),:) + monitorlist(i,:); %!
        if hastags, output.taglist{ubl2num(i)} = [output.taglist{ubl2num(i)} ' ' taglist{i}]; end
    end
    origval = output.valuelist(:,1);  % "real counts"
    output.valuelist(:,1) =          origval     ./ output.monitorlist(:,1) * normval;
    output.valuelist(:,2) = sqrt(max(origval,1)) ./ output.monitorlist(:,1) * normval;  % ** max avoids zero error bars. **
else  % otherwise work with weighted averages
    ind = 1;
    for i=ubl2'         %i is the "bin"-number
        lines = binlist(logical(binlist(:,2)==i),1);  %"lines" are all datapoints that belong to this bin
        [avg, sig] = weightedmean( valuelist(lines,1), valuelist(lines,2) ); %(standard weights 1/sig^2)
        output.coordlist (ind, :) = grid(i,:);
        output.valuelist (ind, :) = [avg, sig];
        if hastags, output.taglist{ind} = [output.taglist{ind} ' ' taglist{lines}]; end
        ind = ind + 1;
    end
end   

if output.polarized
    output.pallist = output.coordlist(:,end);
    output.coordlist = output.coordlist(:,1:(end-1));
end
