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
% 'plotvar', 'varname':     Variable for x-axis in the plot (must be among 'var's). If not given, first one is used.
% 'plot', axhandle :        Give either valid axes handle or 'none' to supress plot. New window, if not given. 
% 'plotstyle', {'s1',..}:   Strings defining the marker style for each pal-series (e.g. 'or', '*b', etc.)
% 'monitor', monval :       Monitor to use for mormalization. If not given, use first M1 value of first scan.
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
% 'constraint', 'p1=2*p2;..': Constraints for the fit parameters. Enter arbitrary number of linear(!) constraints
%                           separated by ";", where 'p1', 'p2' etc. denote the parameters of the fit function.
% Switches:
% 'overplot'    : adds the plot to the current axes.
% 'details'     : show detailed information on the data. 
% 'nobin'       : no binning (default if only one file given)
% 'nooutput'    : suppress all text output except errors. 
% 'nolegend'    : do not put a legend
% 'noplot'      : do not plot the results. (like 'plot none')
% 'showfit'     : Write Fit results in the graphics window
% 'llb'         : Use input routine for LLB scan file format
% '..'          : Use same parameter list as for previous call of nplot
%
% Output: 
%   avgdata:   list of binned and averaged data
%   fitresult: If a fit has ben performed, result of fit parameters

% P. Steffens, 06/2012



%%
% **Not tested for case of PALs and ROIs at the same time

%% Check input
knownoptions = {'var','plotvar','plot','plotstyle','monitor','time','legend','offset','setpal','step','start','end','maxdist','reintegrateimps','only','xtransform','ytransform','calc','fit','startval','fitvar','constraint'};
knownswitches = {'overplot','details','nobin','nooutput','nolegend','noplot','showfit','llb','..'};


% Restore parameters from previous call of nplot if necessary, and store new ones
persistent lastnplotparams
if any(strcmpi(varargin,'..')) && ~isempty(lastnplotparams), 
    lastparamind = find(strcmpi(varargin,'..')); lastparamind = lastparamind(end);
    for j = (length(varargin):-1:(lastparamind+1)), varargin{length(lastnplotparams) + j-1} = varargin{j}; end
    for j = 1:length(lastnplotparams), varargin{lastparamind + j -1} = lastnplotparams{j}; end
end    
lastnplotparams = varargin;

% Test if varagin can be interpreted
[message,unknownopt,multipleopt] = checkoptions(varargin, knownoptions, knownswitches);
if ~isempty(message), fprintf(['Warning(s):\n', message]); end
if ~isempty(multipleopt), fprintf('If you give multiple values for the same options, the last occurence is used.\n'); end
    
if any(strcmpi(varargin,'details')), showdetails=true; else showdetails=false; end  % Detailed output? 
if any(strcmpi(varargin,'nooutput')), nooutput=true; else nooutput=false; end  % Suppress all output?

%% Read all files

if any(strcmpi(varargin,'llb')) % Use LLB Scan Format?
    scans = llbtasread(files,'cells');
    
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
    if isempty(data.variables)
        fprintf('Could not determine scanned variable, plotting agains point number. Please use option "var".\n');
        data.variables = {'PNT'};
    end
else
    if ischar(data.variables), try data.variables = eval(data.variables); catch end; end
    if ~iscell(data.variables), data.variables = {data.variables}; end %ensure cell array
end

%% Eventually adjust file format: Check for special case of "fcu"-counting

fcufile = false;
for ns = 1:length(scans)
    if ~isempty(regexp(upper(scans{ns}.COMND),'\sFCU\s', 'once' ))
    % Change a bit the format of DATA in order to treat it the usual way
    % (one count per line)
        fcufile = true;
        scans{ns}.POLAN = {'fcu Up', 'co', 'fcu Down', 'co'};
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
    elseif fcufile
        fprintf('Error: There seem to be files in fcu-mode and others in standard morde. Cannot mix.\n');
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
        fprintf('Error: cannot find monitor or time values for normalization. Please check!\n Evtl. try to normalize on time instead monitor (use option "time").\n');
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

data.paldeflist = {};
data.polarized = false;
data.multichannel = false; % for multidetectors like IMPS

data.raw = 1;
data.type = 'General scan';
data.coordtype = 'general';
data.expname = [];
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
            fprintf('Error: Trying to combine multidetector data with normal data. I don''t know how to do this.\n'); 
            if nargout, avgdata = []; else clear avgdata; end; return; 
        end
        data.multichannel = true;
        channelname = 'ROI';  % Name of the column that designates the channel number
    elseif isfield(scan,'MULTI') %Flatcone scan
        if (scannr>1) && ~data.multichannel
            fprintf('Error: Trying to combine multidetector data with normal data. I don''t know how to do this.\n'); 
            if nargout, avgdata = []; else clear avgdata; end; return; 
        end
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
            
    elseif data.multichannel
        fprintf('Error: Trying to combine multidetector data with normal data. I don''t know how to do this.\n'); 
        if nargout, avgdata = []; else clear avgdata; end; return;
    end
    

    % Analyze "POLAN"-section (pal file)
    if ~isfield(scan,'POLAN') && ~isempty(data.paldeflist) 
        assignpal = readinput('setpal',varargin,'last');
        scan.DATA.PAL = ones(size(scan.DATA.PNT));
        if isempty(assignpal)
            fprintf('Error: File %d does not contain Polarization info. Use "setpal" option to combine with the others.\n',scannr);
            if nargout, avgdata = []; else clear avgdata; end; return;
        end        
    elseif isfield(scan,'POLAN')
        if (scannr>1) && ~data.polarized, fprintf('Error: Trying to combine non-polarized with polarized data.\n'); if nargout, avgdata = []; else clear avgdata; end; return; end
        data.polarized = true;
        % Analyze the information in POLAN and create (append) the list of
        % PAL-Definitions (paldeflist)
        [data.paldeflist, assignpal] = analyzepal(scan, data.paldeflist);
    end
    
    % Append to lists
    coords = [];
    try
        for ii = 1:length(data.variables)
            if ~isfield(scan.DATA, data.variables{ii}) && ~isfield(scan.DATA, upper(data.variables{ii}))
                fprintf('Error: Could not find variable %s in file %d. Check file format and spelling (incl. upper/lower case).\n', data.variables{ii}, scannr);
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
        data.monitorlist = [data.monitorlist; scan.DATA.(moncolumn)];
        data.valuelist = [data.valuelist; monval * [scan.DATA.CNTS ./ scan.DATA.(moncolumn), sqrt(scan.DATA.CNTS) ./ scan.DATA.(moncolumn)]];
        if data.polarized,      data.pallist    =  [data.pallist; assignpal(scan.DATA.PAL)];    end
        if data.multichannel,   data.channellist = [data.channellist; scan.DATA.(channelname)]; end
        for i=length(data.taglist)+(1:size(coords,1)), data.taglist{i} = scan.FILE; end
    catch
        fprintf('Error: Could not combine file %d with the others. (Check file format!)\n', scannr); 
        if nargout, avgdata = []; else clear avgdata; end; return; 
    end
    
end %Scan loop

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
        data.taglist = {data.taglist{goodlines}};
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
        if any(strcmp(varname,secvar)) && i>1
            % treat special case of qk, ql, en  (as part of dqh)
            stind = find(strcmp(varname,{'QK','QL','EN'}))+1;
            varname = 'QH';
        end
        [st,en] = regexp(cmd, ['(?<=\s+D' varname ')(\s+\-?(\d*\.?\d+|\d+\.?\d*))+']);
        if ~isempty(st) 
            stepadd = str2num(cmd(st(1):en(1))) / 2; %#ok<ST2NM>
        elseif isfield(scan1.STEPS,varname)
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
    fprintf('Error: Please give step sizes for all variables. Use option ''step''.\n'); if nargout, avgdata = []; else clear avgdata; end; return; 
end

if showdetails
    fprintf('Scan variables: '); for i=1:length(data.variables), fprintf('%s ', data.variables{i}); end
    fprintf(['.  The step size used for binning is: ' num2str(gridstep(:)','%g ') '\n']);
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
if any(~good) && ~nooutput, fprintf('** %d data points that are not on the scan path have been discarded.\n', sum(~good)); end


% Retain those points that are in the range and on the scan line
data.coordlist = data.coordlist(inrange&good, :);
data.monitorlist = data.monitorlist(inrange&good, :);
data.valuelist = data.valuelist(inrange&good, :);
if data.polarized, data.pallist = data.pallist(inrange&good); end
if data.multichannel, data.channellist = data.channellist(inrange&good); end
data.taglist = {data.taglist{inrange&good}};
lambdai = lambdai(inrange&good); %(use the lambdai below)



% Check if monitor values are reasonable
% (does not concern normalization on time
if ~strcmpi(moncolumn,'TIME')
    if any(data.monitorlist < data.valuelist(:,1)./data.valuelist(:,2) * 5) && ~nooutput
        fprintf('Warning: Normalization on low monitor values! Normalization on TIME may be more accurate. (Use option "time".)\n');
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
        newgridpoints = [newgridpoints; c*ones(size(gridpoints,1),1), p*ones(size(gridpoints,1),1), gridpoints];
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
    fprintf('Warning: combination of data not successful (cmbavg) ! Go on...\n');
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



% 
% if data.polarized    
%     % for simplicity, use pal as additional coordinate
%     data.coordlist = [data.pallist,data.coordlist]; 
%     palgridpoints=zeros(0,size(gridpoints,2)+1);
%     for p = unique(data.pallist)'
%         palgridpoints = [palgridpoints; p*ones(size(gridpoints,1),1), gridpoints];
%     end    
%     %avgdata = cmbavg(data, 'standard', 'grid', [1,gridstep], 'monitor', monval, 'bin', 0);
%     avgdata = cmbavg(data, 'explicit', palgridpoints, 'monitor', monval, 'bin', [1,inf(size(gridstep))]);
%     avgdata.pallist   = avgdata.coordlist(:,1);
%     avgdata.coordlist = avgdata.coordlist(:,2:end);
%     avgdata.paldeflist = data.paldeflist;
%     avgdata.polarized = true;
% else
%     %avgdata = cmbavg(data, 'standard', 'grid', gridstep, 'monitor', monval, 'bin', 0);
%     avgdata = cmbavg(data, 'explicit', gridpoints, 'monitor', monval, 'bin', inf(size(gridstep)));
%     avgdata.polarized = false;
% end

%% Calculate pal-Combination?
calcstring = upper(readinput('calc',varargin,'last'));
if ~isempty(calcstring) && avgdata.polarized
    
    [avgdata, errorstate] = calclinearcombination(avgdata, calcstring, 'output',~nooutput, 'dist', gridstep);
    if errorstate, if nargout, avgdata = []; else clear avgdata; end; return; end
    
%     [st,en] = regexp(calcstring, 'PAL\d+'); %find all PALx entries in calcstring
%     for i=numel(st):-1:1
%         % Convert calcstring in matlab readable string (PALx --> PAL(:,x))
%         calcstring = [ calcstring(1:st(i)+2) '(:,' calcstring(st(i)+3:en(i)) ')' calcstring(en(i)+1:end) ];
%     end
%     % Find coefficients for each PAL-column
%     PAL = eye(max(avgdata.pallist)); %#ok<NASGU>
%     pcoeff = [];
%     try
%         eval(['pcoeff = ' calcstring ';']);
%         if ~nooutput, fprintf(['** Calculating pal-combination. (Works only for linear combinations.) Coefficients are ' num2str(pcoeff(:)','%g  ') '\n']); end
%     catch
%         em = lasterror;
%         fprintf(['Error during calculation: ''' em.message '''  --> Check input! \n']); if nargout, avgdata = []; else clear avgdata; end; return; 
%     end
%     % Use 'scaledata' and 'subtractdata' to do the linear combination
%     cind = find(pcoeff);
%     comb = scaledata( extractsubset(avgdata,'pallist',cind(1)), pcoeff(cind(1)) );
%     if ~nooutput, fprintf('%s \n', comb.legend); end
%     for i = 2:numel(cind)
%         nextdat = extractsubset(avgdata,'pallist',cind(i));
%         comb = subtractdata( comb, scaledata(nextdat , -pcoeff(cind(i)) ), 'nearest', abs(gridstep) );
%         if ~nooutput, fprintf('%s \n', nextdat.legend); end
%     end
%     avgdata = comb;
%     avgdata.polarized = false;
%     if isfield(avgdata,'taglist'), avgdata = rmfield(avgdata,'taglist'); end
%     avgdata.legend = readinput('calc',varargin,'last');

end

%% Calculate x or y-transform

ytransform = readinput('ytransform',varargin,'last');
if ~isempty(ytransform)
    try
        ytransform(strfind(ytransform, 'y')) = 'Y';
        Y = avgdata.valuelist(:,1);
        dY = avgdata.valuelist(:,2);
        smallval = 1e-4;
        Y = avgdata.valuelist(:,1) + smallval;
        avgdata.valuelist(:,2) = eval(ytransform);
        Y = avgdata.valuelist(:,1) - smallval;
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
        X = avgdata.coordlist;
        avgdata.coordlist = eval(xtransform);
    catch
        fprintf('Error on transformation of x-coordinate. Check statement!\n'); if nargout, avgdata = []; else clear avgdata; end; return; 
    end
end




%% Plot

pstyle = readinput('plotstyle',varargin,'last');
legtext = readinput('legend',varargin,'last');
yoffset = readinput('offset',varargin,'last'); if isempty(yoffset), yoffset = 0; end


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

% If some plotstyle or legends explicitly given, assign them:
if ischar(pstyle), try pstyle = eval(pstyle); catch end; end
if ischar(legtext), try legtext = eval(legtext); catch end; end
for np=1:numplots
   if ~isempty(pstyle),  if iscell(pstyle),  dplot{np}.plotstyle = pstyle{mod(np-1, length(pstyle)) +1}; else dplot{np}.plotstyle=pstyle; end; end
   if ~isempty(legtext), if iscell(legtext), dplot{np}.legend   = legtext{mod(np-1, length(legtext))+1}; else dplot{np}.legend=legtext;   end; end
end

    
if ~isempty(readinput('plotvar',varargin,'last'))
    plotvar = find(strcmpi(data.variables,readinput('plotvar',varargin,'last')));
end
if isempty(plotvar), plotvar = 1; end
    
axhandle = readinput('plot',varargin,'last');
if any(strcmpi(varargin,'overplot')), axhandle = gca; end
if any(strcmpi(varargin,'noplot')), axhandle = 'none'; end

if ~strcmpi(axhandle, 'none')
    
    plot1d(dplot,plotvar,axhandle,varargin{:});

    xlabel(data.variables{plotvar});
    ylabel(['Counts (' moncolumn ' ' num2str(monval) ')']);
    
    title(files);
    if length(data.variables)==4 && all(strcmpi(data.variables,{'QH','QK','QL','EN'})) && size(avgdata.coordlist,1)>0
        if all(abs(avgdata.coordlist(1,1:3)-avgdata.coordlist(end,1:3))<=2*maxdist(1:3)), title({files, ['Q = (', num2str(mean(avgdata.coordlist(:,1:3),1),'%g, %g, %g'), ')']});  %#ok<ALIGN>
        elseif abs(avgdata.coordlist(1,4)-avgdata.coordlist(end,4))<=2*maxdist(4), title({files, ['E = ' num2str(mean(avgdata.coordlist(:,4)),'%g') ]}); 
        end
    end
end

%% Fitting

funcname = readinput('fit',varargin,'last');
if ~isempty(funcname)
    funcname = strtrim(funcname);
    % Read start parameters
    startval = readinput('startval',varargin,'last');
    
    % Test if short form for Gaussian Function
    ngauss = [];
    [st,en] = regexp(upper(funcname),'^GAUSS\d*$');  
    if ~isempty(st)
       ngauss = str2double(funcname(st+5:en)); % number of gaussians
       if isempty(ngauss), ngauss=1; end
       funcname = 'fsum(@const';
       for i=1:ngauss, funcname = [funcname, ',@gaussA']; end
       funcname = [funcname, ')'];
    end
    
    % Simply Linear function?
    islinear = ~isempty(regexp(upper(funcname),'^LINEAR$','once')); 
    
    if funcname(1) ~= '@' && isempty(strfind(funcname,'fsum')) && isempty(strfind(funcname,'fmult')), funcname = ['@',funcname]; end  % Ensure "@" for function handles
    
    % Obtain function information
    eval(['ff=' funcname ';']);
    try
        [val fitresult.paramnames paramnum] = ff([],[]);
        if isempty(startval) && isempty(ngauss) && ~islinear % no startvalues given, and no Gauss function
            fprintf('Please give start parameters for fit! Use option ''startval''.\n'); 
            fprintf('Parameter list: '); for c=1:length(fitresult.paramnames), fprintf('%s ',fitresult.paramnames{c}); end, fprintf('\n'); if nargout, avgdata = []; else clear avgdata; end; return;
        end
    catch
        e=lasterror; msg = e.message;
        fprintf('Error while trying to call function "%s" (error message below). Do not perform fitting.\n%s\n', funcname(2:end), msg);
        funcname=[];
    end    
    
end
if ~isempty(funcname) % continue only if successful until here
    
    % Determine which variables are to be fitted
    varindex = readinput('fitvar',varargin,'last');
    varindex(numel(varindex)+1 : paramnum) = 1; %#ok<NASGU>
    
    % Do the fit for each dataset in dplot
    for np = 1:numplots
        if isempty(dplot{np}), continue; end
        if isempty(startval)     % guess parameters if not given (Gauss and linear only)
            if islinear, startval = startvallinear(dplot{np}.coordlist(:,plotvar), dplot{np}.valuelist(:,1)); end
            if ~isempty(ngauss), startval = startvalgauss(dplot{np}.coordlist(:,plotvar), dplot{np}.valuelist(:,1), ngauss); end
            if isempty(startval), fprintf('Automatic guess of start parameters failed!\n'); if nargout, avgdata = []; else clear avgdata; end; return; end
        end
        % Fitting:       
        if ~nooutput, fprintf('** Fitting dataset %d to function: ', np); end
        [fitval,fitresult.optparam(np,1:paramnum),fitresult.errors(np,1:paramnum),paramoutput] = ...
            funcfit(eval(funcname), dplot{np}.coordlist(:,plotvar), dplot{np}.valuelist(:,1), dplot{np}.valuelist(:,2), startval, varindex,varargin{:}) ; 
%            eval(['funcfit( dplot{np},plotvar,' funcname ',startval,varindex);']);
            
        
        % Plotting:
        if ~strcmpi(axhandle, 'none')
            
            gdata = guidata(gcf);
            
            % Try to find color for line
            if isfield(dplot{np},'plotstyle')
                if isfield(dplot{np}.plotstyle,'color'), linecolor = dplot{np}.plotstyle.color; else linecolor = dplot{np}.plotstyle(end); end
            elseif isfield(gdata,'datahandles') && numel(gdata.datahandles)>=numplots
                linecolor = get(gdata.datahandles(end-numplots+np),'color');
            else
                linecolor = 'k';
            end

            % Calculate fit line and plot into graph   
            xvals = linspace(min(dplot{np}.coordlist(:,plotvar)), max(dplot{np}.coordlist(:,plotvar)), 1000);
            fitline = ff( fitresult.optparam(np,1:paramnum), xvals);     
            ph = plot(xvals, fitline + yoffset, '-');
            set(ph, 'Color', linecolor, 'linewidth', 1, 'tag', 'fitline');
            if isfield(gdata,'fitlinehandles'), gdata.fitlinehandles = [gdata.fitlinehandles, ph]; end
            fitresult.line.x = xvals;
            fitresult.line.y = fitline;
            
            guidata(gcf,gdata);
            
            % Write parameters in the graph window
            if any(strcmpi(varargin,'showfit'))
                xl = xlim(gca); yl = ylim(gca);
                text(xl(1),yl(2),paramoutput,'verticalalignment','top','fontname','courier','tag','fitresults');
            end
        end
    end    
    
end

%%  
plotresult = avgdata;  % Store result in global variable "plotresult", which has to be declared in Matlab-workspace before ('global plotresult')
if nargout==0,  clear avgdata; end

