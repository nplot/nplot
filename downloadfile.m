function success = downloadfile(name,server,directory,username)

% Downloads data file directly from ILL spectrometer computer.
% (Is normally not called directly; use transfer.m or checkfile.m instead.)
% Uses scp on local machine.
% To make it work, the command "scp" must be accessible from the current
% directory (test with "system scp" from Matlab command line), which is normally 
% the case under Linux and Mac; on Windows you may install openssh or
% another package and make sure it is included in the Path variable of the
% system. (Same for "ssh" which is needed by other .m files.)
%
% Instructions to generate a key pair (avoids password for scp, ssh, etc.):
% If you do not already have one (!), type on the local machine:
%       ssh-keygen -t rsa
% (This creates a file .ssh/id_rsa.pub under your home directory)
% Copy this to the remote machine and append to the list of authorized keys:
%       scp .ssh/id_rsa.pub in8@in8:~
%       ssh in8@in8 "cat id_rsa.pub >> .ssh/authorized_keys ; rm id_rsa.pub" 
% Now the connection should be possible without password authorization
% (test by "ssh in8@in8")

% P. Steffens, 09/2013



% [knownservers, directory, defaultserver] = getoption('knownservers','defaultdirectory','defaultserver');

if ~iscell(name), h{1}=name; name=h; clear h; end  % ensure cell array
if nargin < 4, username = server; end
if nargin < 3, directory = '/users/data'; end

filestr = [];
for fn=1:length(name)
    if isempty(directory), filestr = [filestr, name{fn},' '];
    else filestr = [filestr, directory,'/',name{fn},' ']; end
end

%% Do remote copy
% This may depend on the OS. Under Windows, other command may be necessary
% if no scp installed.
[re,ou] = system(['scp ' username '@' server ':"' filestr '" .']);

% in case of error "Filename does not match request": see for ex. https://www.nas.nasa.gov/hecc/support/kb/troubleshooting-scp-file-transfer-failure-with-protocol-error_577.html
% Solution: switch to sftp, or do "scp -T ..."


if strfind(ou,'denied')
    fprintf('Access denied. Retry in new console window (enter password). Exit console manually.\nAvoid this complication by generating a key pair. Type help downloadfile for help.\n');
    [re,ou] = system(['scp ' username '@' server ':"' filestr '" . &']);  %retry same in new window (to allow for password typing). Note: in this case, [re,ou] always [0,''] ...
end

if re ~= 0
    if ~isempty(ou), fprintf('%s\n',ou); end
    if nargout > 0, success = zeros(1,length(name)); end
else
    fprintf('File(s) copied from %s:\n',server);
    j=0;
    for i=1:length(name)
        fprintf('%s ',name{i}); j=j+numel(name{i})+1; if j>80, fprintf('\n'); j=0; end
    end
    fprintf('\n');    
    
    if nargout>0
        % Check for success (check only if filenames present)
        for fn=1:length(name)
        % check if file is there
        [fpath, fname, fext] = fileparts(name{fn}); % split (in case name{fn} contains path)
        success(fn) = ~isempty(dir([fname,fext]));
        end
    end
end


