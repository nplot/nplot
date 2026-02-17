function transfer(varargin)

% function transfer(varargin)
%
% copy files from Serdon to local directory (overwrites local files)
%
% transfer 012345:end    copy all from 012345 to latest
% transfer all           copy all files in data directory
% transfer -n            copy latest n files
% transfer exp 4-01-2345 copy all data from exp 4-01-2345
% transfer [...] from in20
% transfer [...] only ascii/nxs
%
% This is the new version of tranfer.m (replaces old version for direct instrument access)

% P. Steffens 02/2026


instrument = readinput('from',varargin);    if isempty(instrument), instrument = getoption('defaultinstrument'); end
filext = readinput('only',varargin);        if ~isempty(filext) && ~any(strcmpi(filext,{'ascii','nxs'})), error('Allowed values for "only": "ascii" or "nxs"'); end
expnum = readinput('exp',varargin,'noconvert');
allfiles = strcmpi(varargin{1},'all'); % Copy all files?
firstnumtoend = str2double(cell2mat(regexpmatch(varargin{1},'\d+(?=:end)')));   % use of ':end'  (NaN if not given)
newestnum = str2double(cell2mat(regexpmatch(varargin{1},'(?<=^-)\d+$')));       % use of '-n'  

if isempty(expnum)
    if ~allfiles && isnan(newestnum) 
        if isnan(firstnumtoend)
            filelist = multifilename(varargin{1});
            [~, filepath] = finddata(filelist{1},instrument);  % find location of first file
        else
            [~, filepath] = finddata(num2str(firstnumtoend,'%06d'),instrument);  % find location of first file
        end
    else
        [~, filepath] = finddata(instrument); % find actual location
    end
else
    d = dir(['//serdon/illdata/*/' instrument '/exp_' expnum '/rawdata']);
    if ~isempty(d)
        filepath = d(1).folder;
    else
        warning('File location not found.');
        return;
    end
end


d = dir(filepath);
% create file list
if ~allfiles && isnan(newestnum) && isnan(firstnumtoend)
    filelist = multifilename(varargin{1});
else
    numors = cellfun(@(t) cell2mat(regexpmatch(t,'\d+')), {d.name}, 'UniformOutput', false); % extract number part of filename
    numors = cell2mat(cellfun(@(t) str2double(t), numors,  'UniformOutput', false)); % transform to array
    if allfiles
        firstnumor = min(numors);
        lastnumor = max(numors);
%         fprintf('Lowest numor: %06d  Highest numor: %06d\n',firstnumor, lastnumor);
    elseif ~isnan(firstnumtoend)
        firstnumor = firstnumtoend;
        lastnumor = max(numors);
%         fprintf('Highest numor: %06d\n',lastnumor);
    elseif ~isnan(newestnum)
        firstnumor = max(numors) - newestnum +1;
        lastnumor = max(numors);
%         fprintf('Numor: %06d...%06d\n',firstnumor, lastnumor);
    end
    % analyze if all numors in range exist
    datarangestr = '';
    prevexists = false;
    existlist = [];
    n = firstnumor;
    while n <= lastnumor+1
         if any(n == numors)
             existlist = [existlist, n]; %#ok<*AGROW> 
             if ~prevexists
                if  ~isempty(datarangestr), datarangestr = [datarangestr, ', ']; end
                datarangestr = [datarangestr, num2str(n,'%06d')]; 
                rangestart = n;
             end
             prevexists = true;
         else
             if prevexists && n-1>rangestart, datarangestr = [datarangestr, ' .. ' num2str(n-1,'%06d')]; end
             prevexists = false;
         end
         n = n+1;
    end
    fprintf('Data range : %s\n', datarangestr) ;
    % create file list
    filelist = cellfun(@(t) num2str(t,'%06d'), num2cell(existlist), 'UniformOutput', false);
    if isempty(filext)
        filelist = [filelist,  cellfun(@(t) [num2str(t,'%06d') '.nxs'], filelist, 'UniformOutput', false)  ];
    elseif any(strcmpi(filext,{'nxs','nexus'}))
        filelist = cellfun(@(t) [num2str(t,'%06d') '.nxs'], filelist, 'UniformOutput', false);
    end
end

% copy if exists
ncop = 0;
for i=1:length(filelist)
    if any(strcmpi(filelist{i},{d.name}))
        try
            copyfile([filepath filesep filelist{i}]);
            ncop = ncop+1;
        catch
            warning([filelist{i} ' not copied.']);
        end
    end
end
switch ncop
    case 0, fprintf('No file ');
    case 1, fprintf('1 file ');
    otherwise, fprintf('%d files ',ncop);
end
fprintf('copied.\n');



% 
% 
% 
% %% Old version below
% 
% 
% function transfer(filenames, server, directory)
% 
% % function transfer(filenames, server, directory)
% %
% % Copy files from server to local directory (overwrites local files).
% % Server is optional, use defaultserver in option.m otherwise
% % filenames = 012345:end    copy all from 012345 to latest
% % filenames = all           copy all files in data directory
% % By default, transfer.m searches in the currently active directory (/users/data)
% % 
% 
% 
% % P. Steffens, 02/2024
% 
% 
% [defaultdirectory, defaultserver,username] = getoption('defaultdirectory','defaultserver','defaultuser');
% 
% 
% if nargin>1, dwnlsrv = server; else dwnlsrv=defaultserver; end
% if nargin>2, dwnldir = directory; else dwnldir=defaultdirectory; end
% 
% 
% allfiles = strcmpi(filenames,'all'); % Copy all files?
% [st,en] = regexpi(filenames,'\d+(?=:end)');  % Test for use of ':end'
% 
% %%
% 
% 
% 
% % username = dwnlsrv;
% % username = 'nomad';
% 
% filelist = [];
% 
% if ~isempty(st) || allfiles   % using "end" or "all"
%     % list content of data directory
%     [ret,ou]=system(['ssh ' username '@' dwnlsrv ' dir ' dwnldir]);
%     if ret ~= 0 
%         if strfind(ou,'denied'), fprintf('Access denied. You probably need to generate a key pair for this functionality. Type help downloadfile for help.\n');
%         else fprintf('Error on contacting remote host for directory listing.\n');
%         end
%         return;
%     end
%     datafiles = regexpmatch(ou,'\d+');
%     if isempty(datafiles), fprintf('No data files found in %s.\n',dwnldir); return; end
%     lastfilenum=0; firstfilenum=inf;
%     for i=1:length(datafiles)
%         filenum = str2num(datafiles{i});
%         lastfilenum = max([filenum, lastfilenum]);
%         firstfilenum = min([filenum, firstfilenum]);
%     end
%     if allfiles, fprintf('Lowest/highest numor found in %s: %d / %d.\n',dwnldir, firstfilenum, lastfilenum);
%     else fprintf('Highest numor found in %s is %d.\n',dwnldir,lastfilenum);
%     end
%     % create filenames cell array
%     if allfiles
%         filelist = datafiles;
%     else
%         filenum = str2double(filenames(st:en));
%         format = ['%0' num2str(en-st+1) 'd'];  % (nb digits filename)
%         i=1;
%         while filenum <= lastfilenum
%             filelist{i} = num2str(filenum,format);
%             filenum = filenum+1; i = i+1;
%         end
%     end
%     
% else
%     filelist = multifilename(filenames);
%     
% end    
% 
% if isempty(filelist), fprintf('No files to copy.\n'); return; end
%     
% 
% % Try to download from server
% downloadfile(filelist,dwnlsrv,dwnldir);
% 
% 
% 
% %% old version
% 
% % if isempty(st)
% %     filelist = multifilename(filenames);
% %     for j= 1:length(filelist)   
% %         % Try to download from server
% %         [dwnl,dwnlsrv] = downloadfile(filelist{j},dwnl,dwnlsrv);
% %     end
% % else  % download all, until first error appears
% %     success = 0;
% %     filenum = str2double(filenames(st:en));
% %     while success == 0
% %         [dwnl,dwnlsrv,success] = downloadfile(num2str(filenum,'%06d'),dwnl,dwnlsrv);
% %         filenum = filenum+1;
% %     end
% % end