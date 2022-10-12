function [avgdata,fitresult] = nplot(files, varargin)

%
% Syntax: [avgdata,fitresult] = nplot(files, varargin)
%
% Examples:     nplot('0123[45:48,50]', 'monitor', 1000, 'var', {'QH','QK','QL'}, 'calc', '2*pal1-(pal2+pal4)' );
%               nplot 0123[45,67] time 10 plotstyle ob plot gca ;
% The parameter list in varagin can contain an arbitrary combination of pairs (parametername, parametervalue):
% Possible options:
% 'var', {'v1','v2',..} :   which variables to use. If not given, take the scanned variables (of first scan). 
%                           This choice determines the coordinates to be considered, the others are ignored.
%                           For Flatcone data, you can also choose 'twotheta'.
% 'xvar', 'varname':        Variable for x-axis in the plot (must be among 'var's). If not given, first one is used.
% 'yvar', 'varname' :       If you want to plot sth else than CNTS, give the column name here. 
%                           You can also plot values of zeros., param., or varia.-sections of the file header
%                           If this option is used, NO normalization and averaging is done.
% 'plotaxes', axhandle :    Give either valid axes handle or 'none' to supress plot. New window, if not given. 
% 'plotstyle', {'s1',..}:   Strings defining the marker style for each pal-series (e.g. 'or', '*b', etc.)
%                           For multiple, can be like "plotstyle {'or','fb'}" or "plotstyle or|fb"
% 'monitor', monval :       Monitor to use for normalization. If not given, use first M1 value of first scan.
% 'time', time (s) :        If given instead of monitor, normalize on time (in seconds)
% 'legend', legtext:        A text (can be cell array) to display as legend
% 'offset', offset:         Shift plot on y-axis by [offset]
% 'setpal', palnum :        Assign scans without POLAN section to this paldef (if mixed pol/unpol. scans).
% 'step', [s1, s2, ..] :    Stepsize for binning. If not given, try to use half stepsize of first scan.
% 'start', [v1, v2, ..] :   Start point
% 'end', ... :              End point
% 'maxdist', [m1, m2, ..]:  Maximum distance of points to scanline (discard points with larger distance). In 
%                           case of more than one coordinate, the scan path is defined either as the connecting
%                           line of start and end point (if both given) or by the step size.
% 'reintegrateimps', rois:  Reintegrate the IMPS multidetector with the new ROI's (rois= 9x4 matrix). 
%                           Needs to access the corresponding ".multi" file. 
% 'only', which :           retain only selected channels or pal-states. (which) like 'pal1', 'roi[2:4,6]', 'chan31' etc.
% 'calc', 'a*pal1+b*pal2..':Calculate linear combination of pal-sets.
% 'xtransform',expr:        Transform x-coordinate by expr, for example expr = '2*x+0.1' (valid Matlab expression)
% 'ytransform',expr:        Transform y-coordinate by expr, for example expr = 'log(y)+1, expr='y.^2+sqrt(y)', etc.
% 'fit', 'funcname':        Fit function 'funcname' to data. Type 'fitfunctions' for a list.
% 'startval', [v1,v2,..]:   Start values for fit parameters. For funcname='gaussX' (where X=number of Gaussians), 
%                           an automatic guess is performed, if startval not given.
% 'fitvar', [0,1,..]:       Give 1's for parameters to be fitted, 0's for parameters to be held constant.
%                           If not given, all are variable.
% 'common', [1,2,..]:       Indices of common parameters (equal for all datasets if multiple fitting)
% 'constraint', 'p1=2*p2;..': Constraints for the fit parameters. Enter arbitrary number of linear(!) constraints
%                           separated by ";", where 'p1', 'p2' etc. denote the parameters of the fit function.
% Switches:
% 'overplot'    : adds the plot to the current axes.
% 'details'     : show detailed information on the data. 
% 'nobin'       : no binning (default if only one file given)
% 'nooutput'    : suppress all text output except errors. 
% 'nolegend'    : do not put a legend
% 'noplot'      : do not plot the results. (like 'plot none')
% 'showfit'     : Write fit results in the graphics window
% 'globalfit'   : Simultaneous fitting of all datasets
% 'llb'         : Use input routine for LLB scan file format
% 'panda'       : Use input routine for Panda scan file format (FRM-2)
% 'nopalanalysis: Ignore content of POLAN, just use pal-values as in file 
% 'FCSumAll'    : For a scan with Flatcone, sum up all channels
% '..'          : Use same parameter list as for previous call of nplot (parameters can be added, '*' to overwrite previous)
%
% Output: 
%   avgdata     : list of binned and averaged data
%   fitresult   : If a fit has been performed, the resulting fit parameters

% P. Steffens, 07/2019



%%
% **Not tested for case of PALs and ROIs at the same time

%% Check input
knownoptions = {'var','xvar','yvar','plotaxes','plotstyle','monitor','time','legend','offset','setpal','step','start','end','maxdist','reintegrateimps','only','xtransform','ytransform','calc','fit','startval','fitvar','common','constraint'};
knownswitches = {'overplot','details','nobin','nooutput','nolegend','noplot','globalfit','showfit','llb','panda','nopalanalysis','fcsumall','..'};

% Set output options
if any(strcmpi(varargin,'details')), showdetails=true; else showdetails=false; end  % Detailed output? 
if any(strcmpi(varargin,'nooutput')), nooutput=true; else nooutput=false; end  % Suppress all output?
warningstring =[];

% Restore parameters from previous call of nplot if necessary, and store new ones
persistent lastnplotparams
if any(strcmpi(varargin,'..')) && ~isempty(lastnplotparams) 
    lastparamind = find(strcmpi(varargin,'..')); lastparamind = lastparamind(end);
    for j = (length(varargin):-1:(lastparamind+1)), varargin{length(lastnplotparams) + j-1} = varargin{j}; end
    for j = 1:length(lastnplotparams), varargin{lastparamind + j -1} = lastnplotparams{j}; end
end    
lastnplotparams = varargin;


% Check for '*' in input string
% (to override params of previous calls without warning)
notwarnmultiple = find(strncmp('*',varargin,1));
for j=notwarnmultiple, if j>1 && isempty(strfind(lower(varargin{j-1}),'plotstyle')), varargin{j} = varargin{j}(2:end); end; end %#ok<STREMP>


% Replace some old option names by new ones
translate = {'plotvar','xvar'; 'plot','plotaxes'};
for tr = 1:size(translate,1)
    [val,rest] = readinput(translate{tr,1}, varargin);
    if isempty(val), continue; elseif ~iscell(val), val = {val}; end
    if isempty(rest), rest = {}; end 
    varargin = rest; 
    for r=1:length(val)
        varargin = [varargin, translate(tr,2), val(r)];
    end
    if ~isempty(val) && showdetails
        fprintf('Option name %s has been replaced by %s\n', translate{tr,1}, translate{tr,2});
    end
end

% Test if varagin can be interpreted
[message,~,multipleopt] = checkoptions(varargin, knownoptions, knownswitches,notwarnmultiple);
if ~isempty(message), warning(message); end
if ~isempty(multipleopt), fprintf('If you give multiple values for the same options, the last occurence is used.\n(Use * in front of option name to suppress this warning.)\n'); end
if any(strcmpi(varargin,'common')) && ~any(strcmpi(varargin,'globalfit')), fprintf('Common-option only used for simultaneous fitting (switch ''globalfit'').\n'); end

%% Read all files

if any(strcmpi(varargin,'llb')) % Use LLB Scan Format?
    scans = llbtasread(files,'cells');
elseif any(strcmpi(varargin,'panda')) % Use Panda Scan Format?
    scans = tasreadpanda(files,'cells');
       
else
    scans = tasread(files,'download','cells');  % ILL
end


if isempty(scans), return; end

scan1 = scans{1};

%% Evtl reintegrate

newrois = readinput('reintegrateimps',varargin,'last');
if ~isempty(newrois)
    for scannr = 1:length(scans)
        scans{scannr} = reintegrateimps(scans{scannr}, newrois);
    end
end


%% Determine the coordinates to be stored
% Either given in varargin, or inferred from scan command
% Give multiple variables as cell array of strings
data.variables = readinput('var',varargin,'last');

if isempty(data.variables)
%     [st,en] = regexp(upper(scan1.COMND),'(?<=(BS|SC)\s+)\w+');           % ** Allow for multiple variables ?! **
    [st,en] = regexp(upper(scan1.COMND),'(?<=\s+D)\w+(?=\s+\-?(\d*\.?\d+|\d+\.?\d*))');           
    % To recognize scanned variables, look at "D.." parts of COMND (the given steps)
    for i=1:numel(st)
        scanvar = upper(scan1.COMND(st(i):en(i)));
        if ~isempty(scanvar)
            data.variables{length(data.variables)+1} = scanvar;
        end
        if strcmp(scanvar,'QH')
            data.variables{length(data.variables)+1} ='QK'; data.variables{length(data.variables)+1} = 'QL'; data.variables{length(data.variables)+1} ='EN';
        end
    end
    % If xvar option given, check if variable in automatically generated var-list
    xvarname = readinput('xvar',varargin,'last');
    if ~isempty(xvarname) && ~any(strcmpi(xvarname,data.variables))
        data.variables = {xvarname};
        if showdetails
            fprintf('Variable list from scan replaced by %s (argument in xvar option).\n',xvarname);
        end
    end
    if isempty(data.variables)
        fprintf('Could not determine scanned variable, plotting agains point number. Please use option "var".\n');
        data.variables = {'PNT'};
    end
else
    if ischar(data.variables), try data.variables = eval(data.variables); catch, end; end
    if ~iscell(data.variables), data.variables = {data.variables}; end %ensure cell array
end

userealtime = false;
if any(strcmpi(data.variables,'realtime')) % special case, calculate this via getvar.m
    if ~nooutput
        fprintf('Using date and time of measurement as x-variable.\n');
        fprintf('Note that values are approximate only (neglect positioning and pauses) and correspond to end of measured point.\n'); 
    end
    for ns = 1:length(scans)
        scans{ns}.DATA.REALTIME = getvar(scans{ns},'realtime');
        scans{ns}.DATA.columnames = [scans{ns}.DATA.columnames, {'REALTIME'}];
    end
    varargin = [varargin, {'nobin'}]; % do no binning in this case
    userealtime = true;
end


%% Determine also if other column than CNTS is to be on the y-axis

yvar = readinput('yvar',varargin,'last');
if isempty(yvar)
    specialyvar = false;
else
    specialyvar = true;
    varargin = [varargin, {'nobin'}];
    if showdetails
        fprintf('Use column %s for y-axis. Errors are set to NaN.\n', yvar);
    end    
    
    % some guesses to correct simplified input if appropriate
    if ~isfield(scan1.DATA,yvar) && isfield(scan1.DATA,upper(yvar)), yvar = upper(yvar); end  % check if upper case may be better
    if ~isfield(scan1.DATA,yvar) && upper(yvar(1))=='Z' && any(isfield(scan1.ZEROS,{yvar(2:end),upper(yvar(2:end))}))
        yvar = ['ZEROS.' yvar(2:end)];
    end
    
    % test if in zeros, params, or varia, and append value as data column
    if length(yvar)>5 && yvar(6)=='.' && any(strcmpi(yvar(1:5),{'ZEROS','PARAM','VARIA'})) 
        try
            yvarsec=upper(yvar(1:5));
            yvar = yvar(7:end);
            if ~isfield(scan1.(yvarsec),yvar) && isfield(scan1.(yvarsec),upper(yvar)), yvar = upper(yvar); end  % check if upper case may be better
            for ns = 1:length(scans)
                scans{ns}.DATA.([yvarsec,'__',yvar])= scans{ns}.(yvarsec).(yvar) * ones(size(scans{ns}.DATA.PNT));
                scans{ns}.DATA.columnames = [scans{ns}.DATA.columnames, {[yvarsec,'__',yvar]} ];
            end
            yvar = [yvarsec,'__',yvar];
        catch
            fprintf('Error on extracting %s from file header.\n',yvar);
            if nargout, avgdata = []; else clear avgdata; end; return;
        end
    end
        
end

%% Eventually adjust file format: Check for special case of "fcu"-counting

fcufile = false;
for ns = 1:length(scans)
    if isfield(scans{ns},'COMND') && ~isempty(regexp(upper(scans{ns}.COMND),'\sFCU\s', 'once' ))
    % Change a bit the format of DATA in order to treat it the usual way
    % (one count per line)
        fcufile = true;
        scans{ns}.POLAN = {'fcu Up', 'co', 'fcu Down', 'co'};
        if ~isfield(scans{ns}.DATA,'PAL') % if there are already PAL's, do not need the following (for IN22 files)
            fields = fieldnames(scans{ns}.DATA);
            cnr = find(strcmpi('columnames',fields));
            for line =1:size(scans{ns}.DATA.PNT,1)
                for fnr = setdiff(1:length(fields), cnr) % Disregard field "columnnames"
                    if isempty(regexp(upper(fields{fnr}),'UP|DOWN', 'once' ))
                        DATANEW.(fields{fnr})(2*line-1:2*line, 1) = scans{ns}.DATA.(fields{fnr})(line);
                    end
                end
                DATANEW.PAL(2*line-1,1) = 1;
                DATANEW.PAL(2*line,  1) = 2;
                DATANEW.M1(2*line-1,1) = scans{ns}.DATA.M_UP(line);
                DATANEW.M1(2*line,  1) = scans{ns}.DATA.M_DOWN(line);
                DATANEW.CNTS(2*line-1,1) = scans{ns}.DATA.DET_UP(line);
                DATANEW.CNTS(2*line,  1) = scans{ns}.DATA.DET_DOWN(line);
                DATANEW.TIME(2*line-1,1) = scans{ns}.DATA.T_UP(line);
                DATANEW.TIME(2*line,  1) = scans{ns}.DATA.T_DOWN(line);            
            end
            scans{ns}.DATA = DATANEW;
        end
    elseif fcufile
        fprintf(2,'Error: There seem to be files in fcu-mode and others in standard morde. Cannot mix.\n');
        if nargout, avgdata = []; else clear avgdata; end; return;
    end
end
if fcufile && showdetails, fprintf('Input files are in fcu-counting mode.\n'); end
scan1 = scans{1}; %(in case something changed)


%% Determine Normalization
moncolumn = 'M1';
monval = readinput('monitor',varargin,'last');
timeval = readinput('time',varargin,'last'); % if normalized on time
if ~isempty(timeval)
    moncolumn = 'TIME';
    monval = timeval;
end
if isempty(monval)
    if isfield(scan1.PARAM,'TI') 
        moncolumn = 'TIME';
        monval = scan1.PARAM.TI;
    elseif isfield(scan1.PARAM,'MN')
        monval = scan1.PARAM.MN;
    elseif isfield(scan1.DATA,'M1')
        monval = scan1.DATA.M1(1);
    else
        fprintf(2,'Error: cannot find monitor or time values for normalization. Please check!\n Evtl. try to normalize on time instead monitor (use option "time").\n');
    end
end

if showdetails 
    fprintf(['Data are normalized to ' moncolumn ' = ' num2str(monval) '.\n']);
    if isempty(readinput('monitor',varargin,'last')) && isempty(readinput('time',varargin,'last'))
        fprintf('(This is the value found in the first file. Use "time" or "monitor" option to change normalization.)\n');
    end
end

%% Initialize output

global plotresult
global nplotfitresult

data.paldeflist = {};
data.polarized = false;
data.multichannel = false; % for multidetectors like IMPS

data.raw = 1;
data.type = 'General scan';
data.coordtype = 'general';
data.expname = '';
if isfield(scan1,'TITLE'), data.expname = scan1.TITLE; end
data.dataname = files;

data.coordlist   = [];
data.valuelist   = [];
data.monitorlist = [];
data.pallist     = [];
data.channellist = [];
data.taglist     = {};

fitresult = [];



%% Loop over scans to collect all data in one structure

for scannr = 1:length(scans)
    
    scan = scans{scannr};
    
    if isfield(scan.DATA,'ROI')
        if (scannr>1) && ~data.multichannel
            fprintf(2,'Error: Trying to combine multidetector data with normal data. I don''t know how to do this.\n'); 
            if nargout, avgdata = []; else clear avgdata; end; return; 
        end
        data.multichannel = true;
        channelname = 'ROI';  % Name of the column that designates the channel number
    elseif isfield(scan,'MULTI') %Flatcone scan
        if (scannr>1) && ~data.multichannel
            fprintf(2,'Error: Trying to combine multidetector data with normal data. I don''t know how to do this.\n'); 
            if nargout, avgdata = []; else clear avgdata; end; return; 
        end
        if any(strcmpi(varargin,'fcsumall'))
            % take sum of all Flatcone channels; treat like normal scan
            scan.DATA.CNTS = sum(scan.MULTI,2);
        else
            data.multichannel = true;
            actchannel = scan.PARAM.CHAN;
            channelname = 'CHAN';  % Name of the column that designates the channel number
            % Convert MULTI-data in column format to treat in the following
            colform = size(scan.DATA.PNT);
            a4val = getvar(scan,'A4');
            datacolnames = scan.DATA.columnames;
            scan.DATA.CHAN = ones(colform) * actchannel;
            scan.DATA.TWOTHETA = a4val + (actchannel-16)*2.5;
            for ch=[1:(actchannel-1),(actchannel+1):size(scan.MULTI,2)]
                scan.DATA.CHAN = [scan.DATA.CHAN ; ones(colform) * ch];
                scan.DATA.TWOTHETA = [scan.DATA.TWOTHETA; a4val + (ch-16)*2.5];
                for col = 1:length(datacolnames)
                    if ~strcmpi(datacolnames{col},'CNTS')
                        scan.DATA.(datacolnames{col}) = [scan.DATA.(datacolnames{col}); scan.DATA.(datacolnames{col})(1:colform(1))];
                    else
                        scan.DATA.CNTS = [scan.DATA.CNTS; scan.MULTI(:,ch)];
                    end
                end
            end
        end
            
    elseif data.multichannel
        fprintf(2,'Error: Trying to combine multidetector data with normal data. I don''t know how to do this.\n'); 
        if nargout, avgdata = []; else clear avgdata; end; return;
    end
    

    % Analyze "POLAN"-section (pal file)
    if ~isfield(scan,'POLAN') && ~isempty(data.paldeflist) 
        assignpal = readinput('setpal',varargin,'last');
        scan.DATA.PAL = ones(size(scan.DATA.PNT));
        if isempty(assignpal)
            fprintf(2,'Error: File %s does not contain polarization info. Use "setpal" option to combine with the others.\n',scan.FILE);
            if nargout, avgdata = []; else clear avgdata; end; return;
        end        
    elseif isfield(scan,'POLAN')
        if isfield(scan.DATA,'PAL')
            if (scannr>1) && ~data.polarized, fprintf('Error (in %s): Trying to combine non-polarized with polarized data.\n',scan.FILE); if nargout, avgdata = []; else clear avgdata; end; return; end
            data.polarized = true;
            if ~any(strcmpi(varargin,'nopalanalysis'))
                % Analyze the information in POLAN and create (append) the list of
                % PAL-Definitions (paldeflist)
                [data.paldeflist, assignpal] = analyzepal(scan, data.paldeflist);
            else
                assignpal = (1:max(scan.DATA.PAL))'; 
                for ii = (length(data.paldeflist)+1):max(scan.DATA.PAL), data.paldeflist{ii}.PAL = ii; end
            end
        else
            warning('Inconsistent file format in %s. Found POLAN, but no PAL''s. Treat as unpolarized (please check).',scan.FILE);
        end
    end
    
    % Append to lists
    coords = [];
    try
        for ii = 1:length(data.variables)
            if ~isfield(scan.DATA, data.variables{ii}) && ~isfield(scan.DATA, upper(data.variables{ii}))
                fprintf(2,'Error: Could not find variable %s in file %s. Check file format and spelling (incl. upper/lower case).\n', data.variables{ii}, scan.FILE);
                if nargout, avgdata = []; else clear avgdata; end; return; 
            end
            try
                coords = [coords, scan.DATA.(data.variables{ii})];
            catch
                coords = [coords, scan.DATA.(upper(data.variables{ii}))];
            end
        end
        data.coordlist = [data.coordlist; coords];
        if any(scan.DATA.(moncolumn)==0)
            fprintf('Zeros were detected in column %s of file %s, which is used for normalization. You may try normalizing on time by using the ''time'' option.\n',moncolumn,scan.FILE);
            if nargout, avgdata = []; else clear avgdata; end; return;
        end
        if ~isfield(scan.DATA,moncolumn)
            fprintf(2,'Error: Could not find %s (used for normalization) in file %s. Check file format and spelling (incl. upper/lower case).\n', moncolumn, scan.FILE);
            if nargout, avgdata = []; else clear avgdata; end; return; 
        end
        data.monitorlist = [data.monitorlist; scan.DATA.(moncolumn)];
        if ~specialyvar %(use CNTS)
            data.valuelist = [data.valuelist; monval * [scan.DATA.CNTS ./ scan.DATA.(moncolumn), sqrt(scan.DATA.CNTS) ./ scan.DATA.(moncolumn)]];
        else % use other column for y-axis, and NaN's as error
            if ~isfield(scan.DATA,yvar)
                fprintf(2,'Error: Could not find variable %s in file %d. Check file format and spelling (incl. upper/lower case).\n', yvar, scannr);
                if nargout, avgdata = []; else clear avgdata; end; return; 
            end
            if any(strcmpi(yvar,{'M1','M2'}))
                if ~nooutput, warningstring = ['Using normalized ', yvar,' on y-axis and sqrt as error. (check if using divider)']; end
                data.valuelist = [data.valuelist; monval * [scan.DATA.(yvar) ./ scan.DATA.(moncolumn), sqrt(scan.DATA.(yvar)) ./ scan.DATA.(moncolumn)]];
            else
                data.valuelist = [data.valuelist; scan.DATA.(yvar), nan(size(scan.DATA.(yvar)))];
            end
        end
        if data.polarized,      data.pallist    =  [data.pallist; assignpal(scan.DATA.PAL)];    end
        if data.multichannel,   data.channellist = [data.channellist; scan.DATA.(channelname)]; end
        len=length(data.taglist);
        for i=(1:size(coords,1))
            if ~userealtime, data.taglist{len+i} = scan.FILE; 
            else data.taglist{len+i} = [scan.FILE,' (', datestr(scan.DATA.REALTIME(i)), ')'];
            end
        end
    catch
        if scannr==1, fprintf(2,'Error: Problem with file %s. (Check the file format!)\n', scan.FILE); 
        else          fprintf(2,'Error: Could not combine file %s with the others. (There may be a problem with the file format, please check.)\n', scan.FILE); end
        if nargout, avgdata = []; else clear avgdata; end; return; 
    end
    
end %Scan loop

if ~isempty(warningstring), warning(warningstring); warningstring =[]; end %#ok<NASGU>

if showdetails
    fprintf('Input data is ');
    if data.polarized, fprintf('polarized '); else fprintf('not polarized '); end
    if data.multichannel, fprintf('and multi-detector (sort by: %s). ',channelname); else fprintf('and single-detector. '); end
    fprintf('Total number of data points found: %d in %d scans. \n', size(data.coordlist,1), scannr);
end


%% Bin automatically?
if scannr==1
    % if only one scan, then do not bin
    nobinning = true; 
else 
    % otherwise yes, unless switch 'nobin' given
    nobinning = any(strcmpi(varargin,'nobin'));
end


%% If only single channels are selected, do not retain the others

selection = readinput('only', varargin,'last');
if ~isempty(selection)
    try
        [st,en] = regexp(upper(selection),'[A-Z]+');
        selname = upper(selection(st:en));
        eval(['selval = ' selection(en+1:end) ';']);
        if any(strcmpi(selname,{'ROI','CHAN'}))
            goodlines = ismember(data.channellist,selval);
        elseif strcmpi(selname,'PAL')
            goodlines = ismember(data.pallist,selval);
        else 
            fprintf('Bad identifier in ''only'' option. Can not identify column.\n'); 
        end
        data.coordlist = data.coordlist(goodlines, :);
        data.monitorlist = data.monitorlist(goodlines, :);
        data.valuelist = data.valuelist(goodlines, :);
        if data.polarized, data.pallist = data.pallist(goodlines); end
        if data.multichannel, data.channellist = data.channellist(goodlines); end
        data.taglist = data.taglist(goodlines);
        if isfield(data,'dataname'), data.dataname = [selection ': ' data.dataname]; end
        if ~nooutput, fprintf(['Retain only points with %s = ' num2str(selval) '.\n'], selname); end

    catch
        fprintf('Error while evaluation ''only'' option. Ignore it and go on...\n');
    end
end



%% Determine a stepsize

startpoint = readinput('start',varargin,'last');
endpoint = readinput('end',varargin,'last');

gridstep = readinput('step',varargin,'last');

if isempty(gridstep)        % Determine stepsize from scan command (1st scan) 
    for i=1:length(data.variables)
        cmd = upper(scan1.COMND);
        varname = upper(data.variables{i});
        stind = 1;
        secvar = {'QK','QL','EN'}; % "secondary var's"
        if any(strcmp(varname,secvar)) %&& i>1 %**??
            % treat special case of qk, ql, en  (as part of dqh)
            stind = find(strcmp(varname,{'QK','QL','EN'}))+1;
            varname = 'QH';
        end
        [st,en] = regexp(cmd, ['(?<=\s+D' varname ')(\s+\-?(\d*\.?\d+|\d+\.?\d*))+']);
        if ~isempty(st) 
            stepadd = str2num(cmd(st(1):en(1))) / 2; %#ok<ST2NM>
        elseif isfield(scan1,'STEPS') && isfield(scan1.STEPS,varname)
            stepadd = scan1.STEPS.(varname) / 2;
        else
            stind=[]; stepadd=[]; 
        end %stepadd can be array (dqh)
        gridstep = [gridstep, stepadd(stind)]; 
    end
else
    nobinning = false;
    if any(strcmpi(varargin,'nobin')) && ~nooutput, fprintf('Switch ''nobin'' inactive because ''step'' explicitly given.\n'); end
end
plotvar = find(gridstep); % first var with nonzero step becomes plotvar
if ~isempty(plotvar), plotvar = plotvar(1); end
if ~isempty(startpoint) && ~isempty(endpoint)
    gridstep = gridstep(1)/(endpoint(1)-startpoint(1)) * (endpoint - startpoint);
    % Make sure that the stepsize points from start to end
    if showdetails
        fprintf('Start and end point are explicitly given.'); 
        if numel(gridstep)>1, fprintf(' Step size is adapted.'); end 
        fprintf('\n'); 
    end 
end

if isempty(gridstep) && nobinning, gridstep = ones(1,length(data.variables)); end  %(step not used)

if numel(gridstep) ~= length(data.variables) || all(gridstep==0)
    fprintf(2,'Error: Please give step sizes for all variables. Use option ''step''.\n'); if nargout, avgdata = []; else clear avgdata; end; return; 
end

if showdetails
    fprintf('Scan variables: '); for i=1:length(data.variables), fprintf('%s ', data.variables{i}); fprintf('\b'); end
    if nobinning
        fprintf('. No binning of data is performed.\n');
    else
        fprintf(['.  The step size used for binning is: ' num2str(gridstep(:)','%g ') '\n']);
    end
end

%% Discard points that do not belong to the scan
% The data set is N-dim., that means that if N>1 all points do not
% necessarily lie on a line

minlambda = 0;
maxlambda = inf;
if isempty(startpoint), minlambda = -inf; startpoint = data.coordlist(1,:); end
if ~isempty(endpoint), maxlambda = (endpoint-startpoint) * gridstep'; end

maxdist = readinput('maxdist',varargin,'last');
if isempty(maxdist)  % set a limit for the maximum distance of a point to the scan line
    % Problem: the practical precision depends a lot on which variable it is.
    % By default, take 0.01:    
    maxdist = max(gridstep, .01*ones(size(gridstep))); 
    % For variables that are kept explicitly constant, look at the scan to get an idea of the precision:
    if isempty(readinput('step',varargin,'last'))
        for zerostep = find(gridstep==0)
            if isfield(scan1.DATA,data.variables{zerostep})
                maxdist(zerostep) = max(maxdist(zerostep), max(scan1.DATA.(data.variables{zerostep}))-min(scan1.DATA.(data.variables{zerostep})));
            end
        end
    end
end 

if showdetails && numel(maxdist)>1
    fprintf(['The maximum accepted distance to the scan path is ' num2str(maxdist(:)','%g ')]);
    if isempty(readinput('maxdist',varargin,'last')), fprintf('. (Use "maxdist" option to change this.)\n'); else fprintf(' (provided explicitly).\n'); end
end

ndim = length(data.variables);
good = true(size(data.coordlist,1),1);  % index for points that are on the line
inrange = good; % index for points that are between given start and end point

% For each point, determine distance to scanline (defined by startpoint and gridstep)
lambdai = 0;
for nd = 1:ndim
    lambdai = lambdai + (data.coordlist(:,nd) - startpoint(nd)) * gridstep(nd);
end
inrange = inrange & (lambdai > minlambda) & (lambdai < maxlambda);
for nd = 1:ndim
    pproj_nd = startpoint(nd) + lambdai/max(1E-10,sum(gridstep.^2)) * gridstep(nd);
    dist_nd = data.coordlist(:,nd) - pproj_nd;
    good = good & (abs(dist_nd) <= maxdist(nd));
end
if any(~good) && ~nooutput 
    fprintf('** %d data points that are not on the scan path have been discarded.\n', sum(~good)); 
    if showdetails
        fprintf('Scan path defined by (%s) + x*(%s).\n',num2str(startpoint,'%6.4f '),num2str(gridstep,'%6.4f '));
        fprintf('First rejected points are:\n'); disp(num2str(data.coordlist(find(~good,3),:),'%6.4f '));
    end
end


% Retain those points that are in the range and on the scan line
data.coordlist = data.coordlist(inrange&good, :);
data.monitorlist = data.monitorlist(inrange&good, :);
data.valuelist = data.valuelist(inrange&good, :);
if data.polarized, data.pallist = data.pallist(inrange&good); end
if data.multichannel, data.channellist = data.channellist(inrange&good); end
data.taglist = data.taglist(inrange&good);
lambdai = lambdai(inrange&good); %(use the lambdai below)



% Check if monitor values are reasonable
% (does not concern normalization on time
if ~strcmpi(moncolumn,'TIME')
    if any(data.monitorlist < data.valuelist(:,1)./data.valuelist(:,2) * 5) && ~nooutput
        warning('Normalization on low monitor values! Normalization on TIME may be more accurate. (Use option "time".)');
    end
end



%% Binning and  Averaging



% À faire:::
%gridstep (gridstep==0) = .01;
lambdastart = min(lambdai) / (gridstep(:)'*gridstep(:));
lambdaend = max(lambdai) / (gridstep(:)'*gridstep(:));
% The lambdai are the projections of (coord-startpoint) on gridstep


for nd=1:numel(gridstep)  % Construct the points of the "scan" (points to which data are binned)
    gridpoints(:,nd) = (floor(lambdastart):ceil(lambdaend))' * gridstep(nd) + startpoint(nd);
end

% On binning, use channel as additional coordinate. If not present, add it.
% Then do the same for the polarization states
if ~data.multichannel
    data.channellist = zeros(size(data.coordlist,1),1);
end
if ~data.polarized
    data.pallist = zeros(size(data.coordlist,1),1);
end
newgridpoints = zeros(0,size(gridpoints,2)+2);
for p = unique(data.pallist)'
    for c = unique(data.channellist)'
        newgridpoints = [newgridpoints; c*ones(size(gridpoints,1),1), p*ones(size(gridpoints,1),1), gridpoints]; %#ok<*AGROW>
    end
end
data.coordlist = [data.channellist, data.pallist, data.coordlist];
% The first two columns of newgridpoints are now the channel number and the
% pal-state. With this, do a normal averaging:
% (Do not use pal-option of cmbavg)
polarized = data.polarized;
data.polarized = false;

if nobinning
    avgdata = cmbavg(data, 'noavg', 'monitor', monval);
else
    avgdata = cmbavg(data, 'explicit', newgridpoints, 'monitor', monval, 'bin', [1,1,inf(size(gridstep))]);
end

if isempty(avgdata)
    warning('combination of data not successful (cmbavg) ! Go on...');
    avgdata = data;
end

data.polarized = polarized;
% Split again and assign necessary fields of avgdata
if data.multichannel 
    avgdata.channellist = avgdata.coordlist(:,1); 
end
if data.polarized
    avgdata.pallist = avgdata.coordlist(:,2);
    avgdata.paldeflist = data.paldeflist;
end

avgdata.coordlist = avgdata.coordlist(:,3:end);
avgdata.polarized = data.polarized;
avgdata.multichannel = data.multichannel;


%% Calculate pal-Combination?
calcstring = upper(readinput('calc',varargin,'last'));
if ~isempty(calcstring) && avgdata.polarized
    
    [avgdata, errorstate] = calclinearcombination(avgdata, calcstring, 'output',nooutput, 'dist', gridstep);
    if errorstate, if nargout, avgdata = []; else clear avgdata; end; return; end
end

%% Calculate x or y-transform

ytransform = readinput('ytransform',varargin,'last');
if ~isempty(ytransform)
    try
        ytransform(strfind(ytransform, 'y')) = 'Y';
        %Y = avgdata.valuelist(:,1);
        dY = avgdata.valuelist(:,2);
        smallval = 1e-4;
        Y = avgdata.valuelist(:,1) + smallval; %#ok<NASGU>
        avgdata.valuelist(:,2) = eval(ytransform);
        Y = avgdata.valuelist(:,1) - smallval; %#ok<NASGU>
        avgdata.valuelist(:,2) = avgdata.valuelist(:,2) - eval(ytransform);
        avgdata.valuelist(:,2) = avgdata.valuelist(:,2) / (2*smallval) .* dY;
        avgdata.valuelist(:,1) = eval(ytransform);
    catch
        fprintf('Error on transformation of y-coordinate. Check statement!\n'); if nargout, avgdata = []; else clear avgdata; end; return; 
    end        
end

xtransform = readinput('xtransform',varargin,'last');
if ~isempty(xtransform)
    try
        xtransform(strfind(xtransform, 'x')) = 'X';
        X = avgdata.coordlist; %#ok<NASGU>
        avgdata.coordlist = eval(xtransform);
    catch
        fprintf('Error on transformation of x-coordinate. Check statement!\n'); if nargout, avgdata = []; else clear avgdata; end; return; 
    end
end




%% Plot

pstyle = readinput('plotstyle',varargin,'last');
legtext = readinput('legend',varargin,'last');

numplots = 0;  % enumerate the data sets to plot


if avgdata.multichannel 
    whichchan = unique(avgdata.channellist)';
else
    whichchan = nan;
end
    
% Loop over channels and palstates to create individual data sets to plot
for chnum = whichchan
    if isnan(chnum), thisdata = avgdata; else thisdata = extractsubset(avgdata,'channellist',chnum,channelname); end
    if avgdata.polarized
        % split in separate lists, according to pal-state
        for np = unique(thisdata.pallist)'
            numplots = numplots + 1;
            dplot{numplots} = extractsubset(thisdata,'pallist',np);   
        end
    else
        numplots = numplots + 1;
        dplot{numplots} = thisdata;
    end
end

% ggf. interpret pstyle, legtext:
if ischar(pstyle), try pstyle = eval(pstyle); catch, end; end %#ok<*SEPEX>
if ischar(legtext), try legtext = eval(legtext); catch, end; end

% define standard legend text
stdlegtext = files;  posqestring = '';
% for Q- or E-Scans determine what is constant
if length(data.variables)==4 && all(strcmpi(data.variables,{'QH','QK','QL','EN'})) && size(avgdata.coordlist,1)>0
    % use values from scan command rather than read out (positioning error on last digits)
    if isfield(scan1,'COMND')
        [st,en]=regexpi(scan1.COMND,'(?<=qh\s+)([\d\s\.\-]+)(?=dqh)');
        posqe = str2num(scan1.COMND(st:en)); %#ok<ST2NM>
    else posqe = mean(avgdata.coordlist,1);
    end
    if all(max(abs([avgdata.coordlist(:,1)-posqe(1),avgdata.coordlist(:,2)-posqe(2),avgdata.coordlist(:,3)-posqe(3)]),[],1)<=2*maxdist(1:3))
        posqestring = num2str(posqe(1:3),'Q = (%g, %g, %g)'); % Q constant 
    elseif max(abs(avgdata.coordlist(:,4)-posqe(4))) <=2*maxdist(4) % E constant
        posqestring = num2str(posqe(4), 'E = %g');
        if isfield(scan1,'POSQE') && isfield(scan1.POSQE,'UN'), posqestring = [posqestring, ' ', scan1.POSQE.UN]; end
    end
end
if ~isempty(posqestring), stdlegtext = [files, '; ', posqestring]; end

% assign legtext and pstyle to dplot structs
for np=1:numplots
   if ~isempty(pstyle),  if iscell(pstyle),  dplot{np}.plotstyle = pstyle{mod(np-1, length(pstyle)) +1}; else dplot{np}.plotstyle=pstyle; end;    end
   if ~isempty(legtext), if iscell(legtext), dplot{np}.legend   = legtext{mod(np-1, length(legtext))+1}; else dplot{np}.legend=legtext;   end 
   elseif (~isfield(dplot{np},'legend') || isempty(dplot{np}.legend)) && ~any(strcmpi(varargin,'nolegend'))
       dplot{np}.legend  = stdlegtext; % use standard text if no other provided 
   end
end

    
if ~isempty(readinput('xvar',varargin,'last'))
    plotvar = find(strcmpi(data.variables,readinput('xvar',varargin,'last')));
    if isempty(plotvar) 
        fprintf(2,'Error: The variable name provided in "xvar" is not found in the variable list.\n'); 
        if nargout, avgdata = []; else clear avgdata; end; return; 
    end        
elseif isempty(plotvar)
    plotvar = 1; % first by default
end
    
axhandle = readinput('plotaxes',varargin,'last');
if any(strcmpi(varargin,'overplot')), axhandle = gca; end
if any(strcmpi(varargin,'noplot')), axhandle = 'none'; end

if ~strcmpi(axhandle, 'none')
    if isempty(axhandle), figure;     axhandle = axes;  end
    
    for axnum = 1:numel(axhandle) % loop if several axes given
    
    % check if plot is empty
    thesearenewaxes = numel(findobj(get(axhandle(axnum),'children'),'tag','data'))==0; 
    
    % Do plot
    linespecs = plot1d(dplot,'axdim',plotvar,'axhandle',axhandle(axnum),varargin{:});
    
    % Assign linespecs to dplot structs
    for np=1:numplots
%         if ~isfield(dplot{np},'plotstyle')
            dplot{np}.plotstyle = linespecs{np}; 
%         end
    end
    
    % Label x and y axes
    xlabeltext = data.variables{plotvar};
    if ~isempty(xtransform), xlabeltext = strrep(xtransform,'X',data.variables{plotvar}); end
    oldxlabel = get(get(gca,'xlabel'),'string'); % Check if axis already labelled (e.g. from prev. plot)
    if ~isempty(oldxlabel) && ~strcmpi(oldxlabel, xlabeltext) && ~nooutput
        warning('The x-axis label is being changed. Check for consistency.');
    end
    xlabel(xlabeltext);
    
    oldylabel = get(get(gca,'ylabel'),'string'); %(as above for x..)
    if ~specialyvar
        if isempty(ytransform)
            ylabeltext = ['Counts (' moncolumn ' ' num2str(monval) ')'];
        else
            ylabeltext = [strrep(ytransform,'Y','Counts') ' (' moncolumn ' ' num2str(monval) ')'];
        end
    else
        if isempty(ytransform)
            ylabeltext = yvar;
        else
            ylabeltext = strrep(ytransform,'Y',yvar);
        end 
    end
    if ~isempty(oldylabel) && ~strcmpi(oldylabel, ylabeltext) && ~nooutput
        warning('The y-axis label is being changed. Check for consistency (e.g. normalization etc.)');
    end
    ylabel(ylabeltext);
    
    % if using date and time on x-axis, adjust axis labels
    if strcmpi(data.variables{plotvar},'realtime')
        dtvec = datevec(avgdata.coordlist(:,plotvar));
        switch find(max(dtvec,[],1)-min(dtvec,[],1),1,'first') % first value changing (y/m/d/h/min/s)
            case 1
                datetick('x','dd-mmm-yyyy'); xlabel('Date');
            case 2
                if max(dtvec(:,3))-min(dtvec(:,2))>6, datetick('x','dd-mmm-yyyy'); xlabel('Date');
                else datetick('x','ddmmm HH:MM'); xlabel('Date and time'); end
            case 3
                datetick('x','ddmmm HH:MM'); xlabel('Date and time');
            case {4,5,6}
                datetick('x','HH:MM:SS'); xlabel(['Time (on ', datestr(avgdata.coordlist(1,plotvar),'dd mmm yyyy') ,')'] );
        end
    end

    % Set title:
    % (if new window, or if new title compatible with existing one)
    oldtitle = get(get(gca,'title'),'string');
    if  thesearenewaxes || any(strcmpi(oldtitle,files)) 
        newtitle = files; 
    else newtitle=''; 
    end
    if ~isempty(posqestring) && (thesearenewaxes || any(strcmpi(oldtitle,posqestring)))      
        if isempty(newtitle), newtitle = posqestring; else newtitle = {newtitle, posqestring}; end
    end
    
    if ~isempty(newtitle), title(newtitle,'interpreter','none'); else title(''); end
    if thesearenewaxes && isempty(legtext) && numplots<2, legend hide; 
    elseif ~any(strcmpi(varargin,'nolegend')), legend show; end

    end
end

%% Fitting

funcname = readinput('fit',varargin,'last');

if ~isempty(funcname) 
    % Do the fit - either for each dataset in dplot separately, or simultaneously for all
    ind=0;
    for np = 1:numplots
        if isempty(dplot{np}), continue; end
        ind = ind+1;

        flcolor{ind}='r';
        if isfield(dplot{np},'plotstyle')
            if isfield(dplot{np}.plotstyle,'color'), flcolor{ind} = dplot{np}.plotstyle.color;
            elseif any(dplot{np}.plotstyle(end)=='ymcrgbwk'), flcolor{ind} = dplot{np}.plotstyle(end); end
        end
        
        allxdata{ind} = dplot{np}.coordlist(:,plotvar);
        allydata{ind} = dplot{np}.valuelist(:,1); 
        allyerror{ind}= dplot{np}.valuelist(:,2); 
        
        if ~any(strcmpi(varargin,'globalfit'))   % individual fits
        
            if ~nooutput, fprintf('** Fitting dataset %d to function: ', np); end

            fitobj{np} = nfit('graphhandle',axhandle,varargin{:}, 'fitfunction',funcname, ...
                    'xdata', allxdata{ind}, 'ydata', allydata{ind}, 'yerror', allyerror{ind}, ...
                    'linecolor', flcolor{ind});  % (this includes the plot)
                
                
            if isempty(fitobj{np}.optimization)
               warning('No fit performed.');
               continue;
            end
            
            % Fill fitresult output structure
            
            nparam = numel(fitobj{np}.parameters.values);
            fitresult.optparam(np,1:nparam) = fitobj{np}.parameters.values;
            fitresult.errors(np,1:nparam) = fitobj{np}.parameters.errors;
            fitresult.chi2(np) = fitobj{np}.optimization.chi2;

            fitresult.fitvalues = fitobj{np}.fitfunction.call(fitobj{np}.parameters.values,dplot{np}.coordlist(:,plotvar));

            % plus some more, to make it more similar to the (original) nfit objects   ** pass directly nfit object??
            fitresult.parameters = fitobj{np}.parameters;
            if isempty(fitobj{np}.fitline.x)
                fitobj{np}.fitline.x = linspace(min(fitobj{np}.xdata),max(fitobj{np}.xdata),1000);                      % x-values;
                fitobj{np}.fitline.y = fitobj{np}.fitfunction.call(fitobj{np}.parameters.values, fitobj{np}.fitline.x); % y-values
            end
            fitresult.fitline = fitobj{np}.fitline;
        end
        
    end   
    
    if any(strcmpi(varargin,'globalfit'))   % simultaneous fit
        if ~nooutput, fprintf('** Fitting datasets (simultaneously) to function: '); end
        fitobj{1} = nfit('graphhandle', axhandle, varargin{:}, 'fitfunction',funcname, ...
                    'xdata', allxdata, 'ydata', allydata, 'yerror', allyerror, ...
                    'linecolor', flcolor);  % (this includes the plot)
        fitresult.optparam = fitobj{1}.parameters.values;
        fitresult.errors = fitobj{1}.parameters.errors;
        fitresult.chi2 = fitobj{1}.optimization.chi2;
        
        % ** evtl some more (s.o.)
        
    end
    
    
    % Write parameters in the graph window
    if any(strcmpi(varargin,'showfit'))
        xl = xlim(gca); yl = ylim(gca);
        text(xl(1),yl(2),fitobj{end}.optimization.paramoutput,'verticalalignment','top','fontname','courier','tag','fitresults');
    end
    
    % evtl. save nfit objects in figure's guidata
    if ~strcmpi(axhandle, 'none') 
        figdata = guidata(gcf); 
        figdata.nfitobj = fitobj;
        guidata(gcf, figdata);
    end
end

%%  
plotresult = avgdata;  % Store result in global variable "plotresult", which has to be declared in Matlab-workspace before ('global plotresult')
nplotfitresult = fitresult; % Same for fitresult
if nargout==0,  clear avgdata; end

