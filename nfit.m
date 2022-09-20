classdef nfit < handle
    
    % nfit class for data fitting
    %
    % Call:
    % nfit              : Creates new nfit class              
    % nfit(parameters)  : Creates new nfit class with property-value pairs(see below)
    % nfit(axes,...)    : Take data from axes (valid axes handle)
    % nfit(function,...): Use specified fit function
    % nfit(axes,function,...): -"-
    % 
    % Possible parameters:
    % xdata, [xdata]    : x-data (cell array for simultaneous fitting of multiple data)
    % (same for ydata, yerror)
    % graphhandle, ax   : axes used for plotting fit line
    % linecolor, cdef   : color of fitline
    % startval, [vals]  : start parameters for fit
    % fitvar, [1/0]     : array containing "1" for variable and "0" for fixed parameters
    % fitfunction, [func]: function used for fitting
    % background, [func]: function describing baseline
    % constraint, [cstr]: linear constraints (like "p1=2*p2;p9=10;p7+p8=p6" etc.)
    % common, commonind : (for multiple datasets:) indices of common parameters
    % Switches :
    % nodata            : Do not perform automatic search for data in open windows
    % *                 : open interactive window
    % 
    % nfit tries to collect data and errorbars from actual axes, if existent
    % The result of the last call is kept in memory for subsequent use:
    % nfit [methodname] (params): Call nfit method on last nfit object
    % 
    % nfit methods:
    % addfunc(func,varargin)    : add a new function 
    % rmfunc(ncomp)             : remove the function ncomp from the sum of functions
    % addconstraint(constr)     : add (one or several) constraint strings
    % rmconstraint(ncon)        : remove one or several constraints
    % setcommon(commons)        : (for multiple dataset fit:) set indices of common (shared) parameters
    % fix(n)                    : fix one or several parameters
    % clear(n)                  : clear one or several parameters
    % fit                       : perform fit
    % setparam(varargin)        : set parameter values
    % plot                      : plot the function with the current parameters
    % (Type 'help nfit.[methodname]' for help on individual methods)
    %
    % Examples:
    % nfit gauss : Fit a Gaussian to data in currently open window (if possible)
    % nfit gauss * : same, and open interactive window after fitting
    % nfit set p1=1 : set first parameter value
    % nfit fix 1 : fix first parameter
    % nfit plot
    % nfit fit
    % f=nfit(gca): create new object f with data in current axes
    % f=nfit('xdata',xdat,'ydata',ydat,'yerror',err,'fitfunction','linear')
    % f.addfunc(@gaussA); add a gaussian to function definition in f
    % f.setparam([1,1,5,30,.5]) : set parameter values
    % f.plot
    %
    % Type 'fitfunctions' for a list of available functions.
    
    % P. Steffens, 1/2022
    
    properties
        xdata = [];     % Data x-values (n-dim.)
        ydata = [];     % Data y-values (1-d)
        yerror = [];    % Data y errors
        fitfunction = struct('call',[], ...     % function call (funct. handle)
                             'name',[], ...     % funct. name
                             'components',[],...% individual function handles of summands
                             'tag', [] );       % tag for each fct., e.g. BGR
        parameters = struct('values',[],'errors',[],'names',[],'fixed',[],'belongtocomponent',[]);
        constraints = {};
        commonparameters = [];
        graphhandle = [];   % handle to axes containing the data
        fitline    = struct('x',[],'y',[],'handle',[],'plotstyle', struct('color','k','linewidth',2));
        graphelements = {}; % evtl additional elements (bg line, etc.)
%         guihandle = [];     % Handle to gui window
        optimization = [];  % optimization structure
    end
    
    events
        nfitupdate
    end
    
    methods
        
        %% Constructor
        function hobj = nfit(varargin)
            persistent nfitobjmemory;
            %initialize cell arrays
            hobj.fitfunction.components ={};
            hobj.parameters.names= {};
                 
            % Check if working on existing object
            if ~isempty(nfitobjmemory) && nargin>0 && ischar(varargin{1})
                methodnames = methods('nfit');
                shortcut = true;
                switch lower(varargin{1})
                    case '*'        % open interactive window on existing object
                        if length(varargin)==1, nfitgui(nfitobjmemory);
                        else fprintf('Use ''nfit *'' without further parameters.\n'); end %#ok<*SEPEX>
                    case 'fitfunction'          % ** initialize function def
                    case 'background'           % ** inizialize background def                    
                    case 'set'              % set one or several param values (call setparam)
                        nfitobjmemory.setparam(varargin{2:end});
                    case {'fix','clear','rmconstraint'}    
                        % methods calls using readinput to extract one suitable parameter
                        ind = readinput(varargin{1},varargin);
                        nfitobjmemory.(varargin{1})(ind);
                    case methodnames % general method call
                        % first input is a method name of nfit:
                        % call this function with the rest of varargin as arguments
                        nfitobjmemory.(varargin{1})(varargin{2:end});
                    otherwise 
                        shortcut = false;
                end
                if shortcut % one of the above cases has been treated
                    hobj = nfitobjmemory;
                    return; 
                end
            end
            if nargin>0 
                try firstarg = eval(varargin{1});
                catch, firstarg = varargin{1};
                end
            end
            
            
            % Initialize...
            % axes for plotting            
            hobj.graphhandle = readinput('graphhandle',varargin,'last');
            if nargin>0 && numel(firstarg)==1 && ishandle(firstarg) % check if first input is axes handle
                parenthandle = firstarg; % input argument is a handle --> use this to take data from figure
                varargin = varargin(2:end); % Remove first entry from varargin (has been treated)
            else
                fig = get(0,'CurrentFigure');
                if ~isempty(fig), parenthandle = get(fig,'CurrentAxes'); else parenthandle = []; end
            end
           
            % Initialize data :
            % Assign directly, if given
            hobj.xdata = readinput('xdata',varargin,'last');
            hobj.ydata = readinput('ydata',varargin,'last');
            hobj.yerror= readinput('yerror',varargin,'last');
            % check sizes
            if any([iscell(hobj.xdata),iscell(hobj.ydata),iscell(hobj.yerror)])
                if ~all([iscell(hobj.xdata),iscell(hobj.ydata),iscell(hobj.yerror)]) || length(hobj.xdata)~=length(hobj.ydata) || length(hobj.xdata)~=length(hobj.yerror) 
                    fprintf(2,'Error: if cell arrays, xdata, ydata, yerror must have the same length.\n'); return;
                end
                numdata=length(hobj.xdata);
            else numdata =1;
            end
            
            % Otherwise :
            if  isempty(hobj.xdata) && isempty(hobj.ydata) && isempty(hobj.yerror) ...
             && ~any(strcmpi('nodata',varargin)) ... 
             && ~isempty(get(0,'CurrentFigure')) && ~isempty(get(gcf,'CurrentAxes')) % some axes do exist
                figdata = guidata(parenthandle); % try to read guidata
                % If guidata has plotdataset structure (created e.g. by plot1d)
                if isfield(figdata,'plotdataset')
                    try
                        ind = length(figdata.plotdataset);  %** evtl. loop for reading multiple data sets?
                        if isgraphics(parenthandle,'axes')
                            % axes handle provided, find dataset in these axes
                            % If several datasets, take first one found
                            while ind>0 && figdata.plotdataset{ind}.axhandle ~= parenthandle, ind=ind-1; end
                        end
                        hobj.xdata  = figdata.plotdataset{ind}.x;
                        hobj.ydata  = figdata.plotdataset{ind}.y;
                        hobj.yerror = figdata.plotdataset{ind}.dy;
                        if isempty(hobj.graphhandle), hobj.graphhandle = figdata.plotdataset{ind}.axhandle; end
                        hobj.fitline.plotstyle.color = get(figdata.plotdataset{ind}.plothandle,'color');
                    catch, fprintf(2,'Error on reading window''s guidata structure. Exit.\n'); return;
                    end
                else
                    % get xy data directly from graphics object
                    hdat = findobj(parenthandle,'type','hggroup','-or','type','errorbar');
                    if isempty(hobj.graphhandle), hobj.graphhandle = get(hdat,'parent'); end
                    if numel(hdat)>=1 
                        hobj.xdata = get(hdat(1),'xdata');
                        hobj.ydata = get(hdat(1),'ydata');
                        hobj.yerror= get(hdat(1),'udata');
                    end
                end
            end
            
                        
            if iscell(hobj.xdata)
                hobj.fitline = {hobj.fitline};
                for ii=2:length(hobj.xdata), hobj.fitline{ii}=hobj.fitline{1}; end
            end
            
            % evtly. set linecolor
            c = readinput('linecolor',varargin,'last');
            if ~isempty(c) 
                if iscell(c) && iscell(hobj.fitline)
                    for ii=1:length(hobj.fitline), if ii<=length(c),  hobj.fitline{ii}.plotstyle.color = c{ii}; end; end
                elseif iscell(hobj.fitline)
                    for ii=1:length(hobj.fitline), hobj.fitline{ii}.plotstyle.color = c; end
                else
                     hobj.fitline.plotstyle.color = c;
                end
            end

            if isempty(varargin), fprintf('Missing argument on call to nfit.\n'); return; end
            
            % Check if (remaining) first argument is a valid function or abbreviation
            if ischar(varargin{1}) && ~any(strtrim(varargin{1})==' ') && ...
                    (~isempty(makefunction(varargin{1})) || ~isempty(regexpi(varargin{1},'^GAUSS\d*$')) || ~isempty(regexpi(varargin{1},'^LORENTZ\d*$')))
                varargin = [{'fitfunction'}, varargin]; 
            end

            
            funcname = readinput('fitfunction',varargin,'last');
            bgr = readinput('background',varargin,'last');
            if ~isempty(funcname)
                if ischar(funcname), funcname = strtrim(funcname); end
                % Test if short form for certain functions: Gauss, Lorentz
                [st,en] = regexp(upper(funcname),'^GAUSS\d*$');
                if isempty(en), ngauss=0; elseif en==5, ngauss=1; else ngauss = str2double(funcname(st+5:en)); end
                [st,en] = regexp(upper(funcname),'^LORENTZ\d*$'); 
                if isempty(en), nlorentz=0; elseif en==7, nlorentz=1; else nlorentz = str2double(funcname(st+7:en)); end
                % Simply Linear function?
                islinear = ~isempty(regexp(upper(funcname),'^LINEAR$','once'));
                
                if ngauss || nlorentz
                    if isempty(bgr), hobj.addfunc(@const,'as background'); end % else hobj.addfunc(bgr,'as background'); end
                    for n=1:ngauss, hobj.addfunc(@gaussA); end
                    for n=1:nlorentz, hobj.addfunc(@lorentzA); end
                else % Normal function
                    hobj.addfunc(funcname);
                end
                if ~isempty(bgr), hobj.addfunc(bgr,'as background'); end
               
                % Initialize start parameters
                startval = readinput('startval',varargin,'last');
                if ~isempty(startval)
                    if numel(startval) == numdata*size(hobj.parameters.values,2)
                        hobj.parameters.values = reshape(startval',size(hobj.parameters.values,2),numdata)';
                    elseif numel(startval) >= size(hobj.parameters.values,2)
                        hobj.parameters.values = startval(1:size(hobj.parameters.values,2));
                    else fprintf('Number of starting values not sufficient.\n'); 
                    end
                elseif islinear % can guess
                    for nd=1:numdata
                        if ~iscell(hobj.xdata), xd =hobj.xdata; yd = hobj.ydata; else xd=hobj.xdata{nd}; yd=hobj.ydata{nd}; end
                        hobj.parameters.values(nd,:) = startvallinear(xd, yd); 
                    end
                elseif ngauss||nlorentz % can guess
                    for nd=1:numdata
                        if ~iscell(hobj.xdata), xd =hobj.xdata; yd = hobj.ydata; else xd=hobj.xdata{nd}; yd=hobj.ydata{nd}; end
                        hobj.parameters.values(nd,:) = startvalgauss(xd, yd, ngauss+nlorentz);
                    end
                end
                
                % evtl. initialize fixed variables 
                varindex = readinput('fitvar',varargin,'last');
                hobj.parameters.fixed(~varindex) = 1;
            end
            
            % Add evtl. constraints
            hobj.addconstraint(readinput('constraint',varargin));
            
            % Set evtl. common parameters
            hobj.setcommon(readinput('common',varargin));

            
            % If a fit is possible, try to perform one
            if ~isempty(hobj.xdata)  && ~isempty(hobj.ydata) && ~isempty(hobj.yerror) && ~isempty(hobj.fitfunction.call) && all(isfinite(hobj.parameters.values(:)))
                hobj.fit;
            end
            
            nfitobjmemory = hobj;   % save object to persistent variable
%             if nargout==0, clear hobj; end
            
            if any(strcmp('*',varargin))
                nfitgui(hobj);  % If *, open interactive gui window on the new nfit object
            end
               

        end
        
        %% Add a function        
        function hobj = addfunc(hobj, func, varargin)
            % hobj = addfunc(hobj, func, varargin)
            % Add a component to the sum of fit functions
            % func = function handle or string
            % varargin can contain a tag, e.g. 'as background', ...
            [func, msg, paramnames, paramnum, description] = makefunction(func);
            if isempty(func)
                if ~isempty(strfind(msg,'Error'))
                    warning(msg); return;
                elseif ~isempty(msg) 
                    fprintf('%s\n',msg); 
                end
            end
            if iscell(hobj.xdata), datalines=1:length(hobj.xdata); else datalines=1; end
            ind = size(hobj.parameters.values,2) + (1:paramnum);
            hobj.parameters.names(ind) = paramnames;
            hobj.parameters.values(datalines,ind) = nan;
            hobj.parameters.errors(datalines,ind) = nan;
            hobj.parameters.fixed(datalines,ind) = 0;
            hobj.parameters.belongtocomponent(ind) = length(hobj.fitfunction.components)+1;
            if ~isempty(hobj.fitfunction.name), hobj.fitfunction.name = [hobj.fitfunction.name, ' + ']; end
            hobj.fitfunction.name = [hobj.fitfunction.name, description];
            hobj.fitfunction.components{end+1} = func;
            hobj.fitfunction.call = fsum(hobj.fitfunction.components{:});
            hobj.fitfunction.tag{end+1} = readinput('as',varargin);            
            % ** may have to repair constraint string if there are simply numbered constraints (not like p2:2)
            
            hobj.optimization = []; % reset optimization info
            % Notifiy changes to evtl. NfitGui
            notify(hobj,'nfitupdate');
            if nargout==0, clear hobj; end
        end
        
        %% Remove a function        
        function hobj = rmfunc(hobj, ncomp)
            % remove the component ncomp from the sum of functions
            if ischar(ncomp), ncomp=str2num(ncomp); end %#ok<ST2NM>
            oldparamnum=size(hobj.parameters.values,2);
            limit= find(hobj.parameters.belongtocomponent>ncomp,1,'first');
            ind = hobj.parameters.belongtocomponent ~= ncomp; 
            hobj.parameters.values = hobj.parameters.values(:,ind);
            hobj.parameters.errors = hobj.parameters.errors(:,ind);
            hobj.parameters.fixed  = hobj.parameters.fixed(:,ind);
            hobj.parameters.names  = hobj.parameters.names(ind);
            hobj.parameters.belongtocomponent = hobj.parameters.belongtocomponent(ind);
            premind = find(~ind);
            premcount = sum(~ind); % number of params removed
            ind = hobj.parameters.belongtocomponent > ncomp; 
            hobj.parameters.belongtocomponent(ind) = hobj.parameters.belongtocomponent(ind)-1;
            hobj.fitfunction.components = hobj.fitfunction.components([1:(ncomp-1),(ncomp+1):end]);
            hobj.fitfunction.tag = hobj.fitfunction.tag([1:(ncomp-1),(ncomp+1):end]);
            hobj.fitfunction.call = fsum(hobj.fitfunction.components{:});
            hobj.fitfunction.name = [];  % update function name
            for ind = 1:length(hobj.fitfunction.components)
                if ind>1, hobj.fitfunction.name = [hobj.fitfunction.name, ' + ']; end
                [~, ~, ~, description] = hobj.fitfunction.components{ind}([]);
                hobj.fitfunction.name = [hobj.fitfunction.name, description];
            end
            hobj.optimization = []; % reset optimization info
            % repair/remove constraints:
            % The numbers of params after the removed Fct. have to be decreased by premcount
            for ic=length(hobj.constraints):-1:1
                [s,e]=regexp(hobj.constraints{ic},'(?<=\<p\(?)\d+(?=\)?\>)'); % look for numbers in sth like 'p1', 'p(3)' etc
                for is=numel(s):-1:1
                    pnum = str2double(hobj.constraints{ic}(s(is):e(is)));
                    if length(hobj.constraints{ic})<=e(is) || hobj.constraints{ic}(e(is)+1)~=':' % parameter not of form p2:3
                        datidstr = [':', num2str(ceil(pnum/oldparamnum))];
                        pnum = mod(pnum-1, oldparamnum)+1;
                    else datidstr = '';
                    end
                    if  pnum >= limit
                        hobj.constraints{ic} = [hobj.constraints{ic}(1:s(is)-1), num2str(pnum-premcount), datidstr, hobj.constraints{ic}(e(is)+1:end)];
                    elseif ismember(pnum,premind) % this constraint contains a deleted param; remove it entirely
                        hobj.constraints = hobj.constraints([1:ic-1,ic+1:end]);
                        break;
                    end
                end
            end
            % Notifiy changes to evtl. NfitGui
            notify(hobj,'nfitupdate');  
        end
                            
        %% Perform fit
        function hobj = fit(hobj,varargin)
            % Perform fitting, and plot resulting line
            if isempty(hobj.fitfunction.call), fprintf('Fitting impossible: No fit function defined.\n'); if nargout < 1, clear hobj; end, return; end
            % Check if any zero error bars, and avoid
            nocell = ~iscell(hobj.yerror);
            if nocell, hobj.yerror = {hobj.yerror}; end
            zmsg=false;
            for ii=1:length(hobj.yerror)
                zeroerr = hobj.yerror{ii}==0; 
                if any(zeroerr)
                    if any(strcmpi('details',varargin)) && ~zmsg, fprintf('For fitting, zero errorbars are replaced by a finite number (smallest errorbar found in dataset).\n'); zmsg=true; end
                    if all(zeroerr), hobj.yerror{ii}=hobj.yerror{ii}+eps; else hobj.yerror{ii}(zeroerr) = min(hobj.yerror{ii}(~zeroerr)); end
                end
            end
            if nocell, hobj.yerror = hobj.yerror{1}; end
            % constraint string
            constr = []; 
            if ~isempty(hobj.constraints), constr = hobj.constraints{1}; for c=2:length(hobj.constraints), constr=[constr,';',hobj.constraints{c}]; end, end %#ok<AGROW>
            % Fit            
            try
                [~,hobj.parameters.values,hobj.parameters.errors,hobj.optimization.paramoutput,~,hobj.optimization.chi2] = ...
                    funcfit(hobj.fitfunction.call, hobj.xdata, hobj.ydata, hobj.yerror, hobj.parameters.values, ~hobj.parameters.fixed, 'constraint', constr, 'common', hobj.commonparameters, 'fullparamoutput') ; 
            catch err            
                fprintf(2,'Error while fitting: %s (in funcfit.m)-> Exit.\n\n', err.message);
                return;
            end 
            hobj.plot;
            notify(hobj,'nfitupdate');  % Notifiy changes to evtl. NfitGui
            if nargout < 1, clear hobj; end
        end
        
        %% Set parameter value
        function hobj = setparam(hobj,varargin)
            % Set parameter values. Input either:
            % [npar], [vals] : two numeric arrays of same size, indicating indices and values 
            % string(s) like "p2=..."
            % [array] : values of parameters starting from first one
            for ind=1:nargin-1 % Try to transform into numeric (for if strings given)
                if ischar(varargin{ind}), num = str2num(varargin{ind}); if ~isempty(num), varargin{ind} = num; end; end %#ok<ST2NM>
            end  
            if nargin==3 && isnumeric(varargin{1}) && isnumeric(varargin{2})
                hobj.parameters.values(varargin{1}) = varargin{2};
                hobj.parameters.errors(varargin{1}) = nan;
            elseif nargin>1 && isnumeric(varargin{1})
                hobj.parameters.values(1:numel(varargin{1})) = varargin{1}(:)';
                hobj.parameters.errors(1:numel(varargin{1})) = nan;
            else
               for ind=1:nargin-1
                    pstr = strsplit(strtrim(varargin{ind}),'=');
                    if length(pstr)>1, val=str2double(pstr{2}); else val=[]; end 
                    pstr = strsplit(pstr{1}((pstr{1}>='0' & pstr{1}<='9') | pstr{1}==':'), ':');
                    if length(pstr)>1, dnum = str2double(pstr{2}); else dnum=[]; end
                    if dnum, pnum = size(hobj.parameters.values,2)*(dnum-1) + str2double(pstr{1}); else pnum = str2double(pstr{1}); end
                    if isnumeric(pnum) && isnumeric(val)
                        hobj.setparam(pnum, val);
                    else
                        fprintf('Input for setparam: [npar],[vals] or [values] or strings like "p2=...", "p2:3=..."\n');
                    end
                end
            end
            hobj.optimization = []; % reset optimization info          
            notify(hobj,'nfitupdate');  % Notifiy changes to evtl. NfitGui
            hobj.plot;
            % ** if constraints, evtl check if values consistent with it
        end
        
        %% Fix a parameter
        function hobj = fix(hobj,ind)
            % hold one or several parameters fixed
            if isnumeric(ind), hobj.parameters.fixed(intersect(1:end,ind)) = 1; end
            notify(hobj,'nfitupdate');  % Notifiy changes to evtl. NfitGui
        end
        
        %% Unfix a parameter
        function hobj = clear(hobj,ind)
            % release one or several parameters
            if isnumeric(ind), hobj.parameters.fixed(intersect(1:end,ind)) = 0; end
            notify(hobj,'nfitupdate');  % Notifiy changes to evtl. NfitGui
        end
        
        %% Add a constraint
        function hobj = addconstraint(hobj, varargin)
            % Add one or several constraint strings
            for ci = 1:nargin-1
                if ischar(varargin{ci})
                    cc = regexpmatch(varargin{ci},'[^;]+'); % (replaces strsplit(varargin{ci},';');)
                    hobj.constraints = [hobj.constraints, cc];
                end
            end
            notify(hobj,'nfitupdate');  % Notifiy changes to evtl. NfitGui
        end
        
        %% Remove a constraint
        function hobj = rmconstraint(hobj, ncon)
            % Remove one or several constraint strings
            if isnumeric(ncon)
                hobj.constraints = hobj.constraints(setdiff(1:end,ncon));
            end
            notify(hobj,'nfitupdate');  % Notifiy changes to evtl. NfitGui
        end
        
        %% Set common parameters
        function hobj = setcommon(hobj, commonindices)
            if isnumeric(commonindices)
                hobj.commonparameters = commonindices;
            end
            notify(hobj,'nfitupdate');  % Notifiy changes to evtl. NfitGui
        end
        
        %% Plot fit line
        function hobj = plot(hobj, varargin)
            nocell = ~iscell(hobj.xdata);
            if nocell, xdat = {hobj.xdata}; hobj.fitline = {hobj.fitline}; else xdat = hobj.xdata; end
            % Delete evtl. old lines
            for ii=1:length(hobj.fitline), if ishandle(hobj.fitline{ii}.handle), delete(hobj.fitline{ii}.handle); hobj.fitline{ii}.handle = []; end; end
            % Plot result function
            for ndat = 1:length(xdat)
                hobj.fitline{ndat}.x = linspace(min(xdat{ndat}),max(xdat{ndat}),1000);                % x-values;
                hobj.fitline{ndat}.y = hobj.fitfunction.call(hobj.parameters.values(ndat,:), hobj.fitline{ndat}.x); % y-values
                if ishandle(hobj.graphhandle)
                    hold(hobj.graphhandle,'on');
                    % plot line
                    hobj.fitline{ndat}.handle = plot(hobj.graphhandle, hobj.fitline{ndat}.x, hobj.fitline{ndat}.y, 'tag', 'fitline');
                    % set attributes
                    for fn=fieldnames(hobj.fitline{ndat}.plotstyle)'
                        set(hobj.fitline{ndat}.handle, fn{1}, hobj.fitline{ndat}.plotstyle.(fn{1})); 
                    end
                end
            end
            % Plotting of other elements
            for iel=1:length(hobj.graphelements)
                if ishandle(hobj.graphelements{iel}.handle), delete(hobj.graphelements{iel}.handle); end
            end
            hobj.graphelements = {};
            if ishandle(hobj.graphhandle)
                bgcomp = strcmpi('background',hobj.fitfunction.tag);
                if any(bgcomp)
                    bgfunc = fsum(hobj.fitfunction.components{bgcomp});
                    for ndat = 1:length(xdat)
                        bg.x = hobj.fitline{ndat}.x;
                        bg.y = bgfunc(hobj.parameters.values(ndat,ismember(hobj.parameters.belongtocomponent,find(bgcomp))),hobj.fitline{ndat}.x);
                        bg.handle = plot(hobj.graphhandle, bg.x, bg.y, 'tag', 'background');
                        bg.plotstyle = struct('color',hobj.fitline{ndat}.plotstyle.color,'linewidth',1,'linestyle','--');
                        hobj.graphelements{end+1} = bg;
                    end
                end
                    % ** Make a plot of more individual components ... (according "tag")
                % set attributes for all elements
                for iel=1:length(hobj.graphelements)
                    for fn=fieldnames(hobj.graphelements{iel}.plotstyle)'
                        set(hobj.graphelements{iel}.handle, fn{1}, hobj.graphelements{iel}.plotstyle.(fn{1})); 
                    end
                end
            end    
            if nocell, hobj.fitline = hobj.fitline{1}; end
            if nargout < 1, clear hobj; end
        end
        
        %% Show fit info
        function str = display(hobj)
            function adst(s,v)
                str = [str, s];
                if nargin>1, vals = [vals,v]; end
            end
            
            vals={};
            str ='nfit: ';
            if isempty(hobj.fitfunction.components)
                adst('(No Fit Function defined.)\n'); %return;
            else
                adst([hobj.fitfunction.name, ' ']);
                if isempty(hobj.optimization) 
                    adst('(No fit performed.)');
                end
                adst('\n');
                adst('Parameter list:\n');
                for j=1:length(hobj.parameters.names)
                    adst(['(',num2str(hobj.parameters.belongtocomponent(j)),') p',num2str(j),' : ',hobj.parameters.names{j},'\n']);
                end
            end
            if nargout==0, fprintf(str); clear str; end
        end
        
        %% Destructor
        function delete(hobj)
            if ishandle(hobj.graphhandle)
                % ** evtl. delete this nfit from that figure's guidata
            end
        end
    end
    
end
