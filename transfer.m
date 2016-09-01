function transfer(filenames, server)

% function transfer(filenames, server)
%
% Copy files from server to local directory

% P. Steffens, 06/2010


dwnl = [];
if nargin>1, dwnlsrv = server; else dwnlsrv=[]; end


[st,en] = regexpi(filenames,'\d+(?=:end)');  % Test for use of ':end'

if isempty(st)
    filelist = multifilename(filenames);
    for j= 1:length(filelist)   
        % Try to download from server
        [dwnl,dwnlsrv] = downloadfile(filelist{j},dwnl,dwnlsrv);
    end
else  % download all, until first error appears
    success = 0;
    filenum = str2double(filenames(st:en));
    while success == 0
        [dwnl,dwnlsrv,success] = downloadfile(num2str(filenum,'%06d'),dwnl,dwnlsrv);
        filenum = filenum+1;
    end
end