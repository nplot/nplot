function [datafilelist, sortparam] = datafinder(filenames, varargin)

% varagin:  pairs of:
%           sort [var] step [st]         
%           filter [condition]
%  with [var] = varname
%       [condition] like EN=0, TT<5, etc.
% output:   filelist: cell array
%           sortparam: values of sorted variables (e.g. EN values for 'sort EN')

% P. Steffens 10/2014

datafilelist = {}; sortparam = [];

varinput = readinput('sort',varargin);
if ~isempty(varinput) && ~iscell(varinput), sortvar{1}=varinput; else sortvar = varinput; end % ensure cell array

stepinput = readinput('step',varargin);
if ~isempty(stepinput) && ~iscell(stepinput), sortstep{1}=stepinput; else sortstep = stepinput; end % ensure cell array

filter = readinput('filter',varargin);
if iscell(filter) %% if several filters given, make one (with AND)
    filter{1}=['(',filter{1},')'];
    for f=2:length(filter), filter{1}=[filter{1},'&(',filter{f},')']; end 
    filter = filter{1};
end
if ~isempty(filter)
    eq = strfind(filter,'='); % if single '=', make '=='
    for f=numel(eq):-1:1
        if filter(eq(f)+1)~='=' && filter(eq(f)-1)~='=', filter = [filter(1:eq(f)),'=',filter(eq(f)+1:end)]; end
    end
    [st,en] = regexp(filter,'[A-Z]\w*'); % find varnames in 'filter' and replace
    for v=numel(st):-1:1
        filter = [filter(1:st(v)-1),'getvar(scan,''', filter(st(v):en(v)), ''',''includeparam'')',filter(en(v)+1:end)];
    end
end


filelist=multifilename(filenames);     
if isempty(filelist), fprintf('Empty file list\n'); return; end

notfound = {};
fileskip = {};
fileok = false(1,length(filelist));

for fn = 1:length(filelist)
    scan = tasread(filelist{fn});
    fileerror = false;
    for v = 1:length(sortvar)
        val = getvar(scan,sortvar{v},'includeparam'); 
        if isempty(val), notfound = [notfound, sortvar{v}]; fileerror = true; continue; end 
        varvalue(fn,v) = mean(val); varmax(fn,v) = max(val); varmin(fn,v) = min(val); 
    end
    
    if fileerror, fileskip = [fileskip, filelist{fn}]; continue; end 
        
    %check if filter ok and store 
    if isempty(filter) || eval(['all(',filter,')']), fileok(fn)=true; end
end

if ~isempty(notfound) 
    notfound = unique(notfound);
    fprintf('Variables not found in all files: '); fprintf('%s ',notfound{:}); fprintf('\n');
end

if ~isempty(fileskip)
    fileskip = unique(fileskip);
    fprintf('Files skipped: '); fprintf('%s ',fileskip{:}); fprintf('\n');
end

% Check for max-min deviations larger that stepsize, and round means to stepsize
for v=1:length(sortstep)
    fileok(varmax(:,v)-varmin(:,v)>sortstep{v}) = false;
    varvalue(:,v) = round(varvalue(:,v)/sortstep{v}) * sortstep{v}; % round to stepsize
end

% Go on with filtered files
filelist = filelist(fileok);
varvalue = varvalue(fileok,:); % varmax = varmax(fileok,:); varmin = varmin(fileok,:);

% sort and group
[sortparam,~,ic]= unique(varvalue,'rows');

for sortedrow = 1:size(sortparam,1)
    filegroup = sort(filelist(ic'==sortedrow));
    if length(filegroup)>1
        % equal characters in the beginning
        datafilelist{sortedrow} = '';
        cm=1;
        while cm<min([numel(filegroup{1}),numel(filegroup{end})]) && filegroup{1}(cm)==filegroup{end}(cm)
            datafilelist{sortedrow}(end+1)=filegroup{1}(cm);
            cm=cm+1;
        end
        datafilelist{sortedrow}(end+1)='[';
        for f=1:length(filegroup)
            datafilelist{sortedrow} = [datafilelist{sortedrow}, filegroup{f}(cm:end), ','];
        end
        datafilelist{sortedrow}(end)=']';
    else
        datafilelist{sortedrow} = filegroup{1}; %#ok<*AGROW>
    end
        
end




