function [expno, user, local, date, title, command, hasdata] = tasscaninfo(filename)

%Get scan information from TAS data file

expno=[]; user=[]; local=[]; date=[]; title=[]; command=[]; 
hasdata = false;

[FID,message] = fopen(filename,'rt');
if (~isempty(message))
    disp(['An error occured while opening the file ' filename ': ' message ]);
    expno='file error';
    return
end

while 1
    fline=fgetl(FID); %read line
    if fline==-1, break, end
    pos=findstr(fline,'EXPNO');
    if ~isempty(pos), expno=fline(8:end); end
    pos=findstr(fline,'USER');
    if ~isempty(pos), user=fline(8:end); end 
    pos=findstr(fline,'LOCAL');
    if ~isempty(pos), local=fline(8:end); end 
    pos=findstr(fline,'DATE');
    if ~isempty(pos), date=fline(8:end); end 
    pos=findstr(fline,'TITLE');
    if ~isempty(pos), title=fline(8:end); end 
    pos=findstr(fline,'COMND');
    if ~isempty(pos), command=fline(8:end); end
    pos=findstr(fline,'DATA');
    if ~isempty(pos), hasdata = true; end
    if hasdata && ~isempty(command) && ~isempty(title) && ~isempty(date) && ~isempty(local) && ~isempty(user) && ~isempty(expno), break, end     
end

fclose(FID);
       