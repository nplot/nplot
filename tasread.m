function [scans, nscans, nodata]=tasread(filenames, varargin)

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


% P. Steffens 07/12 - 10/2014



filelist = multifilename(filenames);

% Check if files are present (and download them, if varargin contains 'download')
checkfile(filelist,varargin{:});

nscans= length(filelist);
% notloaded = [];
nodata = [];


% Go thorough loop for all files...

nsc = 0;

for j=1:nscans
    
    filename = filelist{j};
   
    scanfile=[];
    cpol = 0;

    [FID,message] = fopen(filename); % 'rt'
    
    if (~isempty(message))
         disp(['An error occured while opening the file ' filename ': ' message ]);
%         notloaded = [notloaded,nsc]; 
        continue;
    end

    filecontent = fread(FID);   % read whole file at once
    fclose(FID);

    endline = [0, find(filecontent==10 | filecontent==13)'];
    endline = setdiff(endline, endline-1); % in case of double line returns, keep only second
    if endline(end)~=numel(filecontent), endline(end+1) = numel(filecontent)+1; end %#ok<AGROW>
    flen = numel(filecontent);
    filecontent = char(filecontent);
    

    firstletter = filecontent(endline(1:end-1)+1);
    for nl = find(firstletter >='A' & firstletter <='Z')'  % look (only) at lines starting with a letter

        fline = filecontent(endline(nl)+1:endline(nl+1)-1)';
        
        if ~isempty(strfind('INSTR:,EXPNO:,USER_:,LOCAL:,FILE_:,DATE_:,TITLE:,COMND:,FORMT:,TYPE_:',fline(1:6)))
            % Take the rest of the line as value. These lines should appear only once.
            section = fline(1:5); section = section(fline(1:5)>='A' & fline(1:5)<='Z'); 
            scanfile.(section) = strtrim(fline(7:end));
            continue;
        end
        
        switch fline(1:6)
            case 'DATA_:' 
                fline = strtrim(filecontent(endline(nl+1)+1:endline(nl+2)-1)'); % read line with col.names
                if isempty(fline), break, end
                scanfile.DATA.columnames = regexpmatch(fline,'\S+');
                data = sscanf(filecontent(endline(nl+2)+1:flen)','%f'); % read whole block (until non-number appears)
                data = reshape(data, length(scanfile.DATA.columnames), [])';
            case 'MULTI:'
                multi = sscanf(filecontent(endline(nl+1)+1:flen)','%d'); % read whole block (until non-number appears)
                scanfile.MULTI = reshape(multi, 31, [])';
            case 'POLAN:'
                cpol = cpol+1;
                scanfile.POLAN{cpol} = strtrim(fline(7:end)); %char(regexpmatch(fline,'(?<=:\s+)\S+.+'));
            otherwise
                if fline(6)==':', section = fline(1:5); section = section(section>='A' & section<='Z');
                else section = regexpmatch(fline,'[A-Z]*(?=_*:)');
                end
                if isempty(section), continue; end
                
                try
                    % go along line, look for '=' and ','
                    stind = 7;
                    cmind = [find(fline(8:end)==',')+7, numel(fline)+1]; 
                    for eqind = find(fline(8:end-1)=='=')+7
                        val = fline(eqind+1: min(cmind(cmind>eqind))-1);                    
                        numval = str2double(val);
                        if ~isfinite(numval), numval= str2num(val); end %#ok<ST2NM>
                        if isempty(numval)
                            scanfile.(section).(strtrim(fline(stind:eqind-1))) = strtrim(val);
                        else
                            scanfile.(section).(strtrim(fline(stind:eqind-1))) = numval;
                        end
                        stind = min(cmind(cmind>eqind))+1;
                    end
                catch
                    fprintf('Parse error in file %s, line: %s\n',filename,fline);
                end
        end
        

        
    end

    nsc = nsc+1;
    
    % properly format the DATA matrix; each column becomes a separate array 
    if isfield(scanfile,'DATA')
        doublecolind = [];
        for jj=1:size(data,2)
            if isfield(scanfile.DATA,scanfile.DATA.columnames{jj}), doublecolind(end+1) = jj; end  %#ok<AGROW>
            scanfile.DATA.(scanfile.DATA.columnames{jj})=data(:,jj); 
        end
        scanfile.DATA.columnames = scanfile.DATA.columnames(setdiff(1:jj,doublecolind));
        
    else % no data in this file
        scanfile.FORMT = [];
        scanfile.DATA  = [];
        if (nsc>1)&&(isfield(scans(1),'MULTI')), scanfile.MULTI = []; end
        nodata = [nodata,nsc]; %#ok<AGROW>
    end
    
    
    % Assign to output list    
    try
        if any(strcmpi(varargin,'cells')) 
            scans{nsc} = scanfile; %#ok<AGROW>
        else
            scans(nsc) = scanfile; %#ok<AGROW>
        end
    catch
%         notloaded = [notloaded, j];
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