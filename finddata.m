function varargout = finddata(varargin)
% Output : [numor, filepath, exp, cycle, filename]

% Search for data on SERDON
%
% input arguments:
% intrumentname (optional)
% filename or numor (optional)
% experiment number (optional)

% P Steffens 02/2026


instrumentname = getoption('defaultinstrument');
searchname = [];
expname = [];

for i=1:nargin 
    if any(strcmpi(varargin{i},{'thales','in8','in20','in12','in22','in3','in14'})) 
        instrumentname = lower(varargin{i}); 
    elseif any(strcmpi(varargin{i},{'internal','internaluse'}))
        expname = 'internalUse';
    elseif any(varargin{i}=='-')
        expname = varargin{i};
        if length(expname)<4 || any(expname(1:4)~='exp_'), expname = ['exp_', expname]; end
    else 
        searchname = varargin{i}; % ** maybe check for filename ??
    end
end




filepath = []; cycle = []; exp = []; filename = []; numor = [];


if isempty(searchname) && isempty(expname) % suche aktuelles Experiment

    d=dir(['//serdon/illdata/data/', instrumentname]);
    maxdnum=[]; lastfile={}; lastfolder = {};
    for ind = find(isfolder(cellfun(@(t) ['//serdon/illdata/data/',instrumentname,'/',t,'/rawdata'], {d.name}, 'UniformOutput', false)))
        dd= dir(['//serdon/illdata/data/',instrumentname,'/', d(ind).name, '/rawdata/??????.nxs']);
        if ~isempty(dd)
            [maxdnum(ind),ii] = max(cell2mat({dd.datenum}));
            lastfile{ind} = dd(ii).name;
            lastfolder{ind} =dd(ii).folder;
        end
    end
    [~,maxind] = max(maxdnum);
    ss = strsplit(d(maxind).name,'exp_');
    exp = ss{end};
    filename = lastfile{maxind};
    filepath = lastfolder{maxind};
    ss=strsplit(filename,'.');
    numor = ss{1};

    if nargout==0
        fprintf('Experiment: %s\nLast numor %s\nData path: %s\n',exp,numor,filepath);
    end

elseif ~isempty(expname)  % Suche Experiment
    
    d = dir(['//serdon/illdata/*/', instrumentname, '/', expname ]);
    if ~isempty(d)
        filepath = d(1).folder;
        dd=dir([filepath filesep 'rawdata' filesep '??????*']);
        [~,firstnumor] = fileparts(dd(1).name);
        [~,lastnumor]  = fileparts(dd(end).name);
        if nargout==0
            fprintf('Data path:  %s\nData range: %s ... %s\n',filepath,firstnumor,lastnumor);
        end
    else
        fprintf('Not found.\n');
    end


elseif ~isempty(searchname)  % Suche File

    if ~any(searchname=='.'), searchname(end+1)='*'; end

    % ordered list of cycles (newest first)
    cyclist ={};
    d=dir('//serdon/illdata');
    for decade = '7890123456'
        cyclist = [cyclist, sort({d(cellfun(@(t) t(1),{d.name})==decade).name})];
    end
    cyclist = cyclist(end:-1:1);
    
    % Look for file
    for cyc = cyclist
        d = dir(['//serdon/illdata/', cyc{1}, '/', instrumentname, '/*/rawdata/', searchname ]);
        if ~isempty(d)
            filepath = d(end).folder;
            date = d(end).date;
            filename = {d.name};
            cycle = cyc{1};
            ss = strsplit(filepath, {'\','/'});
            if any(strcmpi(ss,'rawdata')), exp = ss{find(strcmpi(ss,'rawdata'),1,'last')-1}; ss = strsplit(exp,'exp_'); exp = ss{end}; end
            ss=strsplit(filename{1},'.');
            numor = ss{1};
            break;
        end
    end

    if nargout==0
        if isempty(filepath)
            fprintf('File not found.\n');
        else
            fprintf('File found:\nExperiment  %s\nDate        %s\nData path   %s\n',exp,date,filepath);
        end
    end
end

if nargout 
    varargout = {numor, filepath, exp, cycle, filename};
end