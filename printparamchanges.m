function printparamchanges(paramname, range, varargin)

% example : 
% printparamchanges za1 -100            checks last 100 files
% printparamchanges da -100 local       looks in local folder (otherwise on SERDON)
%


% PS, 02/2026


instrumentname = getoption('defaultinstrument');
if ~isnan(str2double(range)), range = str2double(range); end

if isnumeric(range) && range<0  % last n scans
    filelist = {};
    if any(strcmpi(varargin,'local'))
        d= dir;
        numors = cellfun(@(t) cell2mat(regexpmatch(t,'\d+')), {d.name}, 'UniformOutput', false); % extract number part of filename
        numors = cell2mat(cellfun(@(t) str2double(t), numors,  'UniformOutput', false)); % transform to array
        numors = sort(unique(numors),'descend');
        for n = numors
            if any(strcmp(num2str(n,'%06d.nxs'),{d.name}))
                filelist = [filelist, num2str(n,'%06d.nxs')] ; %#ok<*AGROW> 
            elseif any(strcmp(num2str(n,'%06d'),{d.name}))
                filelist = [filelist, num2str(n,'%06d')] ;
            end
            if length(filelist) >= abs(range), break; end
        end
    else
        % ordered list of cycles (newest first)
        cyclist ={};
        d=dir('//serdon/illdata');
        for decade = '7890123456'
            cyclist = [cyclist, sort({d(cellfun(@(t) t(1),{d.name})==decade).name})];
        end
        cyclist = cyclist(end:-1:1);
        ci = 1; 
        while ci <= length(cyclist) && length(filelist) < abs(range)
            cyc = cyclist{ci};
            d = dir(['//serdon/illdata/', cyc, '/', instrumentname, '/*/rawdata/???*' ]);
            if ~isempty(d)
                numors = cellfun(@(t) cell2mat(regexpmatch(t,'\d+')), {d.name}, 'UniformOutput', false); % extract number part of filename
                numors = cell2mat(cellfun(@(t) str2double(t), numors,  'UniformOutput', false)); % transform to array
                numors = sort(unique(numors),'descend');
                for n = numors
                    if any(strcmp(num2str(n,'%06d.nxs'),{d.name}))
                        filelist = [filelist, [d(strcmp(num2str(n,'%06d.nxs'),{d.name})).folder filesep num2str(n,'%06d.nxs')] ];
                    elseif any(strcmp(num2str(n,'%06d'),{d.name}))
                        filelist = [filelist, [d(strcmp(num2str(n,'%06d'),{d.name})).folder filesep num2str(n,'%06d')] ];
                    end
                    if length(filelist) >= abs(range), break; end
                end
            end
            ci = ci+1;
        end
    end
end

filelist = filelist(end:-1:1);  % start with oldest file
prevval = [];
for filename = filelist
    if strcmp(filename{1}(end-2:end),'nxs')  % open nexus
        file = nxsread(filename{1});
        val = nxsgetvar(file,paramname);
        date = file.start_time;
        if isfield(file,'user') && isfield(file.user,'proposal'), prop = file.user.proposal; else, prop = []; end
    else        % open ascii
        file = tasread(filename{1});
        val = getvar(file,paramname,'includeparam');
        date = file.DATE;
        if isfield(file,'EXPNO'), prop = file.EXPNO; else, prop = []; end
    end
    if isempty(prevval) || abs(prevval - val) > 1e-6
        s = split(filename{1},{'/','\'});
        if (isempty(prop) || prop == ' ') && length(s)>3, prop = s{end-2}; end      % prop from folder name on serdon
        fprintf('%20s %11s %12s   %s =%7.3f\n',date, prop, s{end}, paramname, val);
        prevval = val;
    end
end


