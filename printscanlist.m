function [slist,scantable] = printscanlist(filenames,varargin)

% Print MAD scan info
% Filenames can be entered in Multifilename-Format
% 'printscanlist all' lists all data files in directory.
% Can also use 'printscanlist 012345:end'
%
% possible options:   
%       outfile [filename]      if given, write in text file on disk, otherwise on screen
%       var [varname]           prints values of variable (you can give this option several times).
%                               You may add a rounding precision by adding '~[prec]' to varname. Prints *** if 
%                               variance too large. Use '~mean' to print the mean value instead.
% return value:
%       struct array with scan list
% 
% Example: printscanlist 0123[45:67] var EN var QVERT~.01 var TT~mean outfile scanlist.txt
%          list=printscanlist('0123[45:67','var EN');

% P. Steffens, 02/2014

if nargin<1
    help printscanlist
    return;
end
toend = false;

if strcmpi(filenames,'all')
    files=dir;
    filenames = '[';
    for i=1:length(files)
        if length(files(i).name)==6 && all(files(i).name >= '000000') && all(files(i).name <= '999999')
            filenames = [filenames, files(i).name, ','];
        end
    end
    filenames(end)=']';
    if length(filenames)<2 
        fprintf('No data files found in current directory.\n'); 
        return; 
    end
elseif contains(filenames,':end')
    filenames = strrep(filenames,':end','');
    toend = true;
end

filelist=multifilename(filenames);

if toend
   files=dir;
   lastfile = filelist{end};
   for i=1:length(files)
      if length(files(i).name)==6 && all(files(i).name >= '000000') && all(files(i).name <= '999999') && str2double(files(i).name)>str2double(lastfile)
          filelist = [filelist, files(i).name];
      end
   end
end
     
if isempty(filelist)
    fprintf('Empty file list\n');
    return;
end
   

fid=1;
outfile = readinput('outfile',varargin);
if ~isempty(outfile)
    fid=fopen(outfile,'w');
end

varinput = readinput('var',varargin);
if ~isempty(varinput) && ~iscell(varinput), vars{1}=varinput; else vars = varinput; end % ensure cell array

lastexpno=[]; lastuser=[]; lastlocal=[]; lasttitle=[];

fprintf('\n');
for i=1:length(filelist)
    slist(i).file = filelist{i}; %#ok<*AGROW>
%          if ~isempty(vars) %if need to extract variables, use tasread.
        scan=tasread(filelist{i});
%         scan.VARIA.GU = 0; scan.VARIA.GL = 0; scan.VARIA.A3P = 0;
        if isempty(scan), continue; end
        if isfield(scan,'EXPNO'), expno = scan.EXPNO; else expno = 'none'; end %#ok<*SEPEX>
        if isfield(scan,'LOCAL'), local = scan.LOCAL; else local = 'none'; end
        if isfield(scan,'TITLE'), title = scan.TITLE; else title = 'none'; end
        if isfield(scan,'COMND'), command = scan.COMND; else command = 'none'; end
        if isfield(scan,'DATE'), date = scan.DATE; else date = 'none'; end
        if isfield(scan,'USER'), user = scan.USER; else user = 'none'; end
        hasdata = isfield(scan,'DATA') && ~isempty(scan.DATA);
%     else     %otherwise, simply tasscaninfo (a bit faster)
%         [expno, user, local, date, title, command, hasdata] = tasscaninfo(filelist{i});    
%         if strcmp(expno,'file error'), continue, end
%     end
    slist(i).date = date; slist(i).command = command;
    scantable.date{i}=date; scantable.command{i}=command;
    nlflag=0;
    if ~strcmp(expno,lastexpno),    fprintf(fid,'**Exp.-No.: %s    ',expno);  lastexpno=expno; nlflag=1; end
    if ~strcmp(user,lastuser),      fprintf(fid,'**User: %s     ',user); lastuser=user; nlflag=1; end
    if ~strcmp(local,lastlocal),    fprintf(fid,'**Local: %s     ',local); lastlocal=local;  nlflag=1; end
    if ~strcmp(title,lasttitle),    fprintf(fid,'**Title: %s     ',title); lasttitle=title;  nlflag=1; end
    if nlflag==1, fprintf(fid,'\n'); end
    if i==1
        fprintf(fid,'\n  File    Date               ');
        for v=1:length(vars), fprintf(fid,' %7s',char(vars{v})); end
        fprintf(fid,'   Command');
        fprintf(fid,'\n ------------------------------------------------------------------------\n');
    end
    fprintf(fid,' %7s  %15s ',filelist{i},date);
    if ~hasdata 
        fprintf(fid,' (empty file)\n'); 
        continue; 
    end
    for v=1:length(vars)
        st = regexp(vars{v},'~[\d*\.\d*]'); % Check if a precision is given
        if ~isempty(st)
            prec = str2double(vars{v}((st+1):end));
        else
            prec = .001;
        end
        tild = strfind(vars{v},'~');
        if isempty(tild), thisvar = vars{v}; else thisvar = vars{v}(1:(tild-1)); end
        
        val = getvar(scan,thisvar,'includeparam');
        if isempty(val), continue; end
        if ~ischar(val)
            if strfind(vars{v},'~mean') 
                val = mean(val); 
            else
                %val = round(val*100./max(max(val),1))/100.*max(max(val),1); % round a little bit
                val = round(val/prec)*prec;
                val = unique(val);
            end
            if numel(val)>1
                scantable.(thisvar)(i) = nan;
                val = '***'; 
            else
                scantable.(thisvar)(i) = val;
                val = num2str(val);
            end
        end
        fprintf(fid,' %7s',val);
        slist(i).(thisvar)=val;
    end

    fprintf(fid,'   %s',command);
    
    fprintf(fid,'\n');
end

if fid>1, fclose(fid); end

if nargout==0, clear slist ; end
if nargout > 1
    scantable.file = filelist;
end
