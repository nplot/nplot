function [scans, nscans, nodata]=tasreadpanda(filenames, varargin)

% [scans, nscans, nodata] = tasread(filenames, varargin)
% Flexible load routine for ILL TAS data file. 
% varargin may contain 'cells' (gives cell array output), 
% 'download' (tries to download file from server if not found) 
%
% Can load multiple files simultaneously when filename is given in the
% format like ".../scandirectory/012[300:456]". In this case, an array of
% scans is returned. nscans is the number of scans. nodata are those scans
% without DATA field.
%
% Assumes no particular order or content of the file, but simply creates a struture 
% containing fields for ALL variables and values found in the file. Further
% processing or organization of the data (for instance approriate for polarization 
% analysis, Flatcone, etc.) is thus to be done by calling routines.
%
% 
%
% P. Steffens 07/11
%


function matches = regexpmatch(text,expr)
% Simulates the behavior of regexp(..,..,'match'), which is not implemented
% in older Matlab-Versions
[st,en] = regexp(text,expr);
matches=[];
for ii=1:numel(st)
    matches{ii}=text(st(ii):en(ii));
end
end

%%

filelist = multifilename(filenames);

nscans= length(filelist);
notloaded = [];
nodata = [];
dwnl = [];
dwnlsrv=[];

% Go thorough loop for all files...

nsc = 0;

for j=1:nscans
    
    filename = filelist{j};
   
    scanfile=[];
    indata=0;
    inmulti=0;
    cpol = 0;
    np=[];

    [FID,message] = fopen(filename,'rt');
    
%     if strcmpi(message,'No such file or directory') && any(strcmpi(varargin,'download'))
%         % Try to download from server
%         [dwnl,dwnlsrv] = downloadfile(filename,dwnl,dwnlsrv);
%         [FID,message] = fopen(filename,'rt');
%         if isempty(message), fprintf('%s copied from %s.\n', filename, dwnlsrv); end
%     end     
    
    if (~isempty(message))
         disp(['An error occured while opening the file ' filename ': ' message ]);
        notloaded = [notloaded,nsc];
        continue;
    end

    while 1
        fline=fgetl(FID); %read line
        if fline==-1, break, end

        %Beginning of line
%         if fline(1)>='A' && fline(1)<='Z' %#ok<ALIGN> %if starts with letter
%             section = regexpmatch(fline,'[A-Z]*(?=_*:)');
%         else section=[]; end
        
        if strfind(fline,'### Scan data')
            section = {'DATA'};
        elseif strcmp(fline(1:3),'###'), indata = 0; section=[];
        else section=[]; end

        %if new section, no longer in DATA or MULTI block
        if ~isempty(section), indata=0; inmulti=0; end

        %if in DATA or MULTI block, add line
        if indata==1
            data(dline+1,:)=sscanf(fline(setdiff(1:length(fline),strfind(fline,';'))),'%f',[1,inf]);
            dline=dline+1;
        end


        if isempty(section) && ~isempty(strfind(fline, '  imps   ')), scanfile.impsmode = true; continue; end
        
        if isempty(section), continue; end

        %Begin of DATA or MULTI block
        if strcmp(section,'DATA')
            fline=fgetl(FID); %read line
            if fline==-1, break, end
            scanfile.DATA.columnames = regexpmatch(fline(setdiff(2:length(fline),strfind(fline,';'))),'\S+');
            % replace ':' in columnames
            for ind=1:length(scanfile.DATA.columnames), scanfile.DATA.columnames{ind}(scanfile.DATA.columnames{ind}==':') = '_'; end
            % translate to ILL conventions
            ind = find(strcmpi(scanfile.DATA.columnames,'h'),1,'first'); if ~isempty(ind), scanfile.DATA.columnames{ind} = 'QH'; end
            ind = find(strcmpi(scanfile.DATA.columnames,'k'),1,'first'); if ~isempty(ind), scanfile.DATA.columnames{ind} = 'QK'; end
            ind = find(strcmpi(scanfile.DATA.columnames,'l'),1,'first'); if ~isempty(ind), scanfile.DATA.columnames{ind} = 'QL'; end
            ind = find(strcmpi(scanfile.DATA.columnames,'e'),1,'first'); if ~isempty(ind), scanfile.DATA.columnames{ind} = 'EN'; end
            ind = find(strcmpi(scanfile.DATA.columnames,'timer'),1,'first'); if ~isempty(ind), scanfile.DATA.columnames{ind} = 'TIME'; end
            ind = find(strcmpi(scanfile.DATA.columnames,'mon1'),1,'first'); if ~isempty(ind), scanfile.DATA.columnames{ind} = 'M1'; end
            ind = find(strcmpi(scanfile.DATA.columnames,'mon2'),1,'first'); if ~isempty(ind), scanfile.DATA.columnames{ind} = 'M2'; end
            ind = find(strcmpi(scanfile.DATA.columnames,'det2'),1,'first'); if ~isempty(ind), scanfile.DATA.columnames{ind} = 'CNTS'; end
            fline=fgetl(FID); %read line
            if fline==-1, break, end
            indata=1;
            dline=0;
            if isempty(np) || ~isfinite(np) , np=0; end;
            data=zeros(np,length(scanfile.DATA.columnames)); 
        end
        if strcmp(section,'MULTI')
            inmulti=1;
            mline=0;
            scanfile.MULTI=zeros(np,31); 
        end

        %If line contains variables and values, extract these and create
        %respective fields of output
        names   = regexpmatch(fline,'([A-Z]\w*)(?=\s*=)');
%         numbers = regexpmatch(fline,'(?<==\s*)([\d.-]+|[a-zA-Z\*]+)'); %         Read only one number each
        numbers = regexpmatch(fline,'(?<==\s*)([\d.-]+(\s*[\d.-]+)*|[a-zA-Z\*]+)');  % Allow for multiple numbers (case for ROIs) 
        for i=1:length(names)
            if all(numbers{i}<='9'), numbers{i}=str2num(numbers{i}); end
            if ~isempty(numbers{i}), scanfile.(section{1}).(names{i})=numbers{i}; end
        end

        if ~isempty(strfind('INSTR,EXPNO,USER,LOCAL,FILE,DATE,TITLE,COMND,FORMT',section{1}))
            % Lines that (should) appear only once
            scanfile.(section{1}) = char(regexpmatch(fline,'(?<=:\s+)\S+.+'));
            if strcmp(section{1},'COMND'), np=str2double(regexpmatch(fline,'(?<=NP\s)\d+')); end
        end
        
        if strcmpi('POLAN',section{1})
            % there are usually several POLAN lines
            cpol = cpol+1;
            scanfile.POLAN{cpol} = char(regexpmatch(fline,'(?<=:\s+)\S+.+'));
        end
        
        
        
        scanfile.FILE = filename;
        
    end
    
    fclose(FID);
    
    nsc = nsc+1;
    
    if isfield(scanfile,'MULTI'), scanfile.MULTI=scanfile.MULTI(1:mline,:); end
    
    % properly format the DATA matrix; each column becomes a separate array 
    if isfield(scanfile,'DATA')
        data=data(1:dline,:);
        for jj=1:size(data,2)
            scanfile.DATA.(scanfile.DATA.columnames{jj})=data(:,jj); 
        end
        
    else 
        scanfile.FORMT = [];
        scanfile.DATA  = [];
        if (nsc>1)&&(isfield(scans(1),'MULTI')), scanfile.MULTI = []; end
        nodata = [nodata,nsc];
    end
    
    
    
    startpunkt = [scanfile.DATA.QH(1),scanfile.DATA.QK(1),scanfile.DATA.QL(1),scanfile.DATA.EN(1)];
    schrittweite = ([scanfile.DATA.QH(end),scanfile.DATA.QK(end),scanfile.DATA.QL(end),scanfile.DATA.EN(end)] - startpunkt) / (numel(scanfile.DATA.QH)-1); 
    scanfile.COMND = ['bs qh ', num2str(startpunkt,'%f '), ' dqh ', num2str(schrittweite,'%f ')];

    % create empty fields
    if ~isfield(scanfile,'PARAM'), scanfile.PARAM = []; end
    
    
    % Assign to output list    
    try
        if any(strcmpi(varargin,'cells')) 
            scans{nsc} = scanfile;
        else
            scans(nsc) = scanfile;
        end
    catch
        notloaded = [notloaded, j];
        fprintf(['Did not load ' scanfile.FILE ': format not identical to previous files in list.\n']);
        nodata = setdiff(nodata,nsc);
        nsc = nsc-1;
    end
   
end

nscans = nsc;
if nscans == 0
    scans = [];
end


end %function