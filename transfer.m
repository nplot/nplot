function transfer(filenames, server, directory)

% function transfer(filenames, server, directory)
%
% Copy files from server to local directory (overwrites local files).
% Server is optional, use defaultserver in option.m otherwise
% filenames = 012345:end    copy all from 012345 to latest
% filenames = all           copy all files in data directory
% By default, transfer.m searches in the currently active directory (/users/data)
% 


% P. Steffens, 02/2024


[defaultdirectory, defaultserver,username] = getoption('defaultdirectory','defaultserver','defaultuser');


if nargin>1, dwnlsrv = server; else dwnlsrv=defaultserver; end
if nargin>2, dwnldir = directory; else dwnldir=defaultdirectory; end


allfiles = strcmpi(filenames,'all'); % Copy all files?
[st,en] = regexpi(filenames,'\d+(?=:end)');  % Test for use of ':end'

%%



% username = dwnlsrv;
% username = 'nomad';

filelist = [];

if ~isempty(st) || allfiles   % using "end" or "all"
    % list content of data directory
    [ret,ou]=system(['ssh ' username '@' dwnlsrv ' dir ' dwnldir]);
    if ret ~= 0 
        if strfind(ou,'denied'), fprintf('Access denied. You probably need to generate a key pair for this functionality. Type help downloadfile for help.\n');
        else fprintf('Error on contacting remote host for directory listing.\n');
        end
        return;
    end
    datafiles = regexpmatch(ou,'\d+');
    if isempty(datafiles), fprintf('No data files found in %s.\n',dwnldir); return; end
    lastfilenum=0; firstfilenum=inf;
    for i=1:length(datafiles)
        filenum = str2num(datafiles{i});
        lastfilenum = max([filenum, lastfilenum]);
        firstfilenum = min([filenum, firstfilenum]);
    end
    if allfiles, fprintf('Lowest/highest numor found in %s: %d / %d.\n',dwnldir, firstfilenum, lastfilenum);
    else fprintf('Highest numor found in %s is %d.\n',dwnldir,lastfilenum);
    end
    % create filenames cell array
    if allfiles
        filelist = datafiles;
    else
        filenum = str2double(filenames(st:en));
        format = ['%0' num2str(en-st+1) 'd'];  % (nb digits filename)
        i=1;
        while filenum <= lastfilenum
            filelist{i} = num2str(filenum,format);
            filenum = filenum+1; i = i+1;
        end
    end
    
else
    filelist = multifilename(filenames);
    
end    

if isempty(filelist), fprintf('No files to copy.\n'); return; end
    

% Try to download from server
downloadfile(filelist,dwnlsrv,dwnldir);



%% old version

% if isempty(st)
%     filelist = multifilename(filenames);
%     for j= 1:length(filelist)   
%         % Try to download from server
%         [dwnl,dwnlsrv] = downloadfile(filelist{j},dwnl,dwnlsrv);
%     end
% else  % download all, until first error appears
%     success = 0;
%     filenum = str2double(filenames(st:en));
%     while success == 0
%         [dwnl,dwnlsrv,success] = downloadfile(num2str(filenum,'%06d'),dwnl,dwnlsrv);
%         filenum = filenum+1;
%     end
% end