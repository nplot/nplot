function fileok = checkfile(filelist, varargin)

% Checks if files in filelist are accessible. Possibility to try download, if not.
% varargin: 
% 'download'     : Try automatic download for files which are not found.
% 'server',serv  : Name of computer for file download
% 'directory',dir: Name of directory for file download
% 'username',user: User name...
% 'noconfirm'    : Don't ask before downloading
%
% fileok: true for accessible files (after evtl download)
%
% If files not found in default directory on server, uses "datafile find" on the
% remote machine to locate missing files. On these, also use gunzip to decompress .gz files

% P. Steffens 09/2013 - 02/2024

% ** question about "./" in front of datafile. When is it necessary?


if ~iscell(filelist), h{1}=filelist; filelist=h; clear h; end  % ensure cell array


%% Check if files exist in local directory

fileok = true;

for fn=1:length(filelist)
    % check if file is there
    fileok(fn) = ~isempty(dir(filelist{fn})); %#ok<*AGROW>
    if fileok(fn), continue; end
    
    [pathstr, name, ext] = fileparts(filelist{fn});
    filename0(1:(6-length(name))) = '0';
    tentname = fullfile(pathstr,[filename0,name,ext]);
    if ~isempty(dir(tentname))
        fprintf('File %s not found, but %s is present.\n',filelist{fn},tentname);
        filelist{fn} = tentname;
        fileok(fn) = true;
    end
end

if all(fileok) || ~any(strcmpi(varargin,'download')), return; end
% Return anyway, if no download desired

%% Try download from default data directory

missingfiles = {};
for fn=find(~fileok), missingfiles = [missingfiles(:)', {filelist{fn}}]; end


[defaultdirectory, defaultserver, defaultuser] = getoption('defaultdirectory','defaultserver','defaultuser');

server = readinput('server',varargin);          if isempty(server), server = defaultserver; end
directory = readinput('directory',varargin);    if isempty(directory), directory = defaultdirectory; end
username = readinput('username',varargin);      if isempty(username), if isempty(readinput('server',varargin)), username = defaultuser; else username = server; end, end

if ~any(strcmpi(varargin,'noconfirm'))
    answer = input(['Try to download missing file(s) from server ' server ' into current local directory?\nPress Enter to continue. (''n'' for No, ''c'' for change server)'], 's');
    if strcmpi(answer,'n'), return; 
    elseif strcmpi(answer,'c') 
        server = input('Enter server name: ','s'); 
        username = input(['Enter user name for login on ',server,': '],'s'); 
        fprintf('Note: You can change the default server and username in your options.m file.\n');
    end
end

success = downloadfile(missingfiles,server,directory);

fileok(~fileok) = success;

if all(success)
    return;
elseif any(success) 
    fprintf('Some missing files were downloaded, but not all.\n');
    return;
end

%% Search elsewhere on server

% if download from standard data dir was not successful (nothing found), try
% to look elsewhere. Use datafile script on instrument computer to locate
% the files.

searchstr = [];
for fn=find(~fileok), searchstr = [searchstr, './datafile find ', filelist{fn}, '; ']; end

[re,ou]=system(['ssh ',username,'@',server,' ', searchstr]);

if strfind(ou,'denied')
    fprintf('Was not able to search in other paths on the server (access denied; type help downloadfile for creating a key pair).\n');   
    return
end

if isempty(ou), fprintf('Nothing found.\n'); return; end

% otherwise, go on, and ask if download these files
fprintf('In the current data directory, the requested files have not been found, but\nat the following locations on %s there are files matching your search:\n%s\n',server,ou);

if ~any(strcmpi(input('Do you want to use these files? ([y]/n)','s'),{'y',''})), return; end

% Create file list and check for double filenames
try
    loadlist = regexpmatch(ou,'\S+');
%     loadlist = loadlist(1:(end-1));
catch
    fprintf('Error on evaluating output from "datafile" script.\n'); return;
end

clear pathstr name ext
notload = false(1,length(loadlist));
for fn=1:length(loadlist)
    [pathstr{fn}, name{fn}, ext{fn}] = fileparts(loadlist{fn});
    if any(strcmpi(name(1:(fn-1)),name{fn}))
        fprintf('Warning: Filename conflict for %s (found same file numor several times)\n',name{fn});
        notload(fn) = true;
    else 
        i = find(strcmp(filelist,name{fn})); % index in original list
        if isempty(i) && ~isempty(ext{fn}), i = find(strcmp(filelist,[name{fn},ext{fn}])); end
        if numel(i)==1, fnindex(fn)=i; else fnindex(fn)=nan; end %?
    end
end
if any(notload), fprintf('Using first occurrence. (To use other files, please transfer manually.)\n'); end

% Download missing files
success = downloadfile(loadlist(~notload),server,'');

fileok(fnindex(~notload)) = success;

% if some zipped files, unzip
gzfiles = strcmpi(ext,'.gz');
for fn=find(gzfiles)
    try
        gunzip([name{fn},'.gz']);
        delete([name{fn},'.gz']);
    catch
        fprintf('Could not unzip %s\n',[name{fn},'.gz']);
        fileok(fnindex(fn)) = 0;
    end
end



     




    
    




