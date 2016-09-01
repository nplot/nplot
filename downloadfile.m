function [answerall, server, success] = downloadfile(name,answerall,server)

% Downloads data file directly from ILL spectrometer computer

% P. Steffens, 06/2010



[knownservers, directory, defaultserver] = getoption('knownservers','defaultdirectory','defaultserver');

%%
if nargin<3 || isempty(server), server = defaultserver; end

answer = [];
if nargin<2 || isempty(answerall) || ~answerall
    answer = input(['Try to download file ' name ' from server ' server ' into current local directory?  (y)es, (n)o, (a)ll, (c)hange server [a] '], 's');
end
if isempty(answer), answer = 'a'; end
answer = strtrim(answer);
answerall = strcmpi(answer,'a');
if strcmpi(answer(1),'n') || ~any(strcmpi(answer(1),{'y','a','c'})), return; end

%% Choose server

if strcmpi(answer,'c')
    fprintf('Available servers: ');
    for i=1:length(knownservers), fprintf('(%d) %s   ',i,knownservers{i}); end 
    server = input('. Enter number:','s');
    servnum = str2double(server); 
    if isempty(servnum)  || servnum<1 || servnum>length(knownservers) , fprintf('Server not known!\n'); return; end
    server = knownservers{servnum};
    answerall = true;
    fprintf('Note: you can change the default server in the options file. (Type "edit options".) \n');
end

%% Username: change here if necessary

username = server;

%% Do copy

if strfind(computer,'WIN')
    password = [username username];
    success = system(['winscp3 ' username ':' password '@' server ':' directory ' /command "get ' name '" exit']);
    % Works for Windows if WinSCP3.exe is included in PATH
else
    success = system(['scp ' username '@' server ':' directory '/' name '  .']);
    % No password necessary if key pair exchanged:
    % Generate on local machine by: ssh-keygen -t rsa
    % (should be in /Users/IN20/.ssh/id_rsa.pub (on MAC))
    % Copy to remote machine as (or append to): .ssh/authorized_keys  (in /users/in20/)
end

    

