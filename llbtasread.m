function scans = llbtasread(filenames, varargin)

% Read LLB triple axis file format.
% Attention, this is a very rudimentary script. It is not sufficientkly
% general to work in every case !!!


filelist = multifilename(filenames);

for i=1:length(filelist), filelist{i} = [filelist{i},'.asc']; end  % Add .asc to filenames


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
    np=0;

    [FID,message] = fopen(filename,'rt');
    scanfile.FILE = filename;
    
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

    
    fline=fgetl(FID); %read line
    [st,en] = regexp(fline(2:end),'\w+');
    if fline(st(10)+1)=='m', moncount=true; else moncount=false; end % Count on M or Ti
    if fline(st(11)+2)=='i', scanfile.PARAM.FX = 1; else scanfile.PARAM.FX = 2; end % ki or kf const
        
    fline=fgetl(FID); %read line
    % Here are all the step sizes
    mind = find(fline=='-'); % Add spaces before "-" signs
    for m=numel(mind):-1:1, fline = [fline(1:mind(m)-1), ' ', fline(mind(m):end)]; end    
    scancom = str2num(fline(2:end));
    qstart = scancom(1:4);
    scanstep = scancom(5:8);
    montime = scancom(10);
    scanfile.PARAM.KFIX = scancom(11);
    scanfile.COMND = ['bs qh ', num2str(qstart), ' dqh ', num2str(scanstep), ' np ', num2str(scancom(9))];
    if moncount, scanfile.COMND = [scanfile.COMND, ' mn ', num2str(montime)]; else [scanfile.COMND, ' ti ', num2str(montime)]; end
    
    
    fline = fgetl(FID);
    fline = fgetl(FID);
    
    scanfile.DATA.QH = [];
    scanfile.DATA.QK = [];
    scanfile.DATA.QL = [];
    scanfile.DATA.EN = [];
    scanfile.DATA.CNTS = [];
    
    while 1
         fline = fgetl(FID);
         if fline==-1, break, end
         fline = str2num(fline);
         
         scanfile.DATA.QH = [scanfile.DATA.QH; qstart(1) + np*scanstep(1)];
         scanfile.DATA.QK = [scanfile.DATA.QK; qstart(2) + np*scanstep(2)];
         scanfile.DATA.QL = [scanfile.DATA.QL; qstart(3) + np*scanstep(3)];
         scanfile.DATA.EN = [scanfile.DATA.EN; qstart(4) + np*scanstep(4)];
         scanfile.DATA.CNTS = [scanfile.DATA.CNTS; fline(2)];
            
         np=np+1;
    end
    scanfile.DATA.NP = ones(size(scanfile.DATA.QH));
    if moncount
        scanfile.DATA.M1 = montime * ones(size(scanfile.DATA.QH));
    else
        scanfile.DATA.TIME = montime * ones(size(scanfile.DATA.QH));
    end
    
    fclose(FID);
    nsc = nsc+1;
    
    
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
