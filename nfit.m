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
    % xdata, [xdata]    : x-data
    % (same for ydata, yerror)
    % startval, [vals]  : start parameters for fit
    % fitvar, [1/0]     : array containing "1" for variable and "0" for fixed parameters
    % fitfunction, [func]: function used for fitting
    % background, [func]: function describing baseline
    % constraint, [cstr]: linear constraints (like "p1=2*p2;p9=10;p7+p8=p6" etc.)
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
    % fix(n)                    : fix one or several parameters
    % clear(n)                  : clear one or several parameters
    % rmconstraint(ncon)        : remove one or several constraints
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
    
    % P. Steffens, 2/2015
    
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
                        else fprintf('Use ''nfit *'' without further parameters.\n'); end
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
            
            % Initialize data :
            % Assign directly, if given
            hobj.xdata = readinput('xdata',varargin,'last');
            hobj.ydata = readinput('ydata',varargin,'last');
            hobj.yerror= readinput('yerror',varargin,'last');
            % Otherwise :
            if  isempty(hobj.xdata) && isempty(hobj.ydata) && isempty(hobj.yerror) ...
                 && ~any(strcmpi('nodata',varargin)) ... 
                 && ~isempty(get(0,'CurrentFigure')) && ~isempty(get(gcf,'CurrentAxes')) % some axes do exist
                % check if first input is axes handle
                if nargin>0 && numel(firstarg)==1 && ishandle(firstarg)                     
                    parenthandle = firstarg; % input argument is a handle --> use this to take data from figure
                    varargin = varargin(2:end); % Remove first entry from varargin (has been treated)
                else
                    parenthandle = gca;
                end
                figdata = guidata(parenthandle); % try to read guidata
                % If guidata has plotdataset structure (created e.g. by plot1d)
                if isfield(figdata,'plotdataset')
                    try
                        ind = length(figdata.plotdataset); 
                        if strcmpi(get(parenthandle,'type'),'axes')
                            % axes handle provided, find dataset in these axes
                            while ind>0 && figdata.plotdataset{ind}.axhandle ~= parenthandle, ind=ind-1; end
                        end
                        hobj.xdata  = figdata.plotdataset{ind}.x;
                        hobj.ydata  = figdata.plotdataset{ind}.y;
                        hobj.yerror = figdata.plotdataset{ind}.dy;
                        hobj.graphhandle = figdata.plotdataset{ind}.axhandle;
                        hobj.fitline.plotstyle.color = get(figdata.plotdataset{ind}.plothandle,'color');
                    catch, fprintf('Error on reading window''s guidata structure. Exit.\n'); return;
                    end
                else
                    % get xy data directly from graphics object
                    hdat = findobj(parenthandle,'type','hggroup');
                    hobj.graphhandle = get(hdat,'parent');
                    if numel(hdat)==1 
                        hobj.xdata = get(hdat,'xdata');
                        hobj.ydata = get(hdat,'ydata');
                        hobj.yerror= get(hdat,'udata');
                    end
                end
            end
            
            % Check if (remaining) first argument is a valid function or abbreviation
            if ~isempty(varargin) && (~isempty(makefunction(varargin{1})) || ~isempty(regexpi(varargin{1},'^GAUSS\d*$')) || ~isempty(regexpi(varargin{1},'^LORENTZ\d*$')))
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
                    if isempty(bgr), hobj.addfunc(@const,'as background'); else hobj.addfunc(bgr,'as background'); end
                    for n=1:ngauss, hobj.addfunc(@gaussA); end
                    for n=1:nlorentz, hobj.addfunc(@lorentzA); end
                else % Normal function
                    hobj.addfunc(funcname);
                end
                if ~isempty(bgr), hobj.addfunc(bgr,'as background'); end
               
                % Initialize start parameters
                startval = readinput('startval',varargin,'last');
                if ~isempty(startval)
                    if numel(startval) >= numel(hobj.parameters.values)
                        hobj.parameters.values = startval(1:numel(hobj.parameters.values));
                    else fprintf('Number of starting values not sufficient.\n'); 
                    end
                elseif islinear % can guess
                     hobj.parameters.values = startvallinear(hobj.xdata, hobj.ydata); 
                elseif ngauss||nlorentz % can guess
                    hobj.parameters.values = startvalgauss(hobj.xdata, hobj.ydata, ngauss+nlorentz);
                end
                
                % evtl. initialize fixed variables 
                varindex = readinput('fitvar',varargin,'last');
                hobj.parameters.fixed(~varindex) = 1;
            end
            
            % Add evtl. constraints
            hobj.addconstraint(readinput('constraint',varargin));
            
            % If a fit is possible, try to perform one
            if ~isempty(hobj.xdata)  && ~isempty(hobj.ydata) && ~isempty(hobj.yerror) && ~isempty(hobj.fitfunction.call) && all(isfinite(hobj.parameters.values))
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
                if ~isempty(msg), fprintf('%s\n',msg); end
                return;
            end
            ind = numel(hobj.parameters.values) + (1:paramnum);
            hobj.parameters.names(ind) = paramnames;
            hobj.parameters.values(ind) = nan;
            hobj.parameters.errors(ind) = nan;
            hobj.parameters.fixed(ind) = 0;
            hobj.parameters.belongtocomponent(ind) = length(hobj.fitfunction.components)+1;
            if ~isempty(hobj.fitfunction.name), hobj.fitfunction.name = [hobj.fitfunction.name, ' + ']; end
            hobj.fitfunction.name = [hobj.fitfunction.name, description];
            hobj.fitfunction.components{end+1} = func;
            hobj.fitfunction.call = fsum(hobj.fitfunction.components{:});
            hobj.fitfunction.tag{end+1} = readinput('as',varargin);            
            hobj.optimization = []; % reset optimization info
            % Notifiy changes to evtl. NfitGui
            notify(hobj,'nfitupdate');
            if nargout==0, clear hobj; end
        end
        
        %% Remove a function        
        function hobj = rmfunc(hobj, ncomp)
            % remove the component ncomp from the sum of functions
            if ischar(ncomp), ncomp=str2num(ncomp); end %#ok<ST2NM>
            limit= find(hobj.parameters.belongtocomponent>ncomp,1,'first');
            ind = hobj.parameters.belongtocomponent ~= ncomp; 
            hobj.parameters.values = hobj.parameters.values(ind);
            hobj.parameters.errors = hobj.parameters.errors(ind);
            hobj.parameters.names  = hobj.parameters.names(ind);
            hobj.parameters.fixed  = hobj.parameters.fixed(ind);
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
                    if  pnum >= limit
                        hobj.constraints{ic} = [hobj.constraints{ic}(1:s(is)-1), num2str(pnum-premcount), hobj.constraints{ic}(e(is)+1:end)];
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
            if isempty(hobj.fitfunction.call), fprintf('Fitting impossible: No fit funciton defined.\n'); if nargout < 1, clear hobj; end, return; end
            % Check if any zero error bars, and avoid
            zeroerr = hobj.yerror==0; 
            if any(zeroerr)
                if any(strcmpi('details',varargin)), fprintf('For fitting, zero errorbars are replaced by a finite number (smallest errorbar found in dataset).\n'); end
                if all(zeroerr), hobj.yerror=hobj.yerror+eps; else hobj.yerror(zeroerr) = min(hobj.yerror(~zeroerr)); end
            end
            % Fit            
            try
                [~,hobj.parameters.values,hobj.parameters.errors,hobj.optimization.paramoutput,~,hobj.optimization.chi2] = ...
                    funcfit(hobj.fitfunction.call, hobj.xdata, hobj.ydata, hobj.yerror, hobj.parameters.values, ~hobj.parameters.fixed, 'constraint', strjoin(hobj.constraints,';')) ; 
            catch err            
                fprintf('Error while fitting: %s. -> Exit.\n', err.message);
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
                try
                    for ind=1:nargin-1
                        pstr = strsplit(strtrim(varargin{ind}),'=');
                        hobj.setparam(str2double(pstr{1}(pstr{1}>='0' & pstr{1}<='9')), str2double(pstr{2}));
                    end
                catch
                	fprintf('Input for setparam: [npar],[vals] or [values] or strings like "p2=..."\n');
                end
            end
            hobj.optimization = []; % reset optimization info          
            notify(hobj,'nfitupdate');  % Notifiy changes to evtl. NfitGui
            % ** if constraints, evtl check if values consistent with it
        end
        
        %% Fix a parameter
        function hobj = fix(hobj,ind)
            % hold one or several parameters fixed
            if isnumeric(ind), hobj.parameters.fixed(intersect(1:end,ind)) = 1; end
        end
        
        %% Unfix a parameter
        function hobj = clear(hobj,ind)
            % release one or several parameters
            if isnumeric(ind), hobj.parameters.fixed(intersect(1:end,ind)) = 0; end
        end
        
        %% Add a constraint
        function hobj = addconstraint(hobj, varargin)
            % Add one or several constraint strings
            for ci = 1:nargin-1
                if ischar(varargin{ci})
                    cc = strsplit(varargin{ci},';');
                    hobj.constraints = [hobj.constraints, cc];
                end
            end
        end
        
        %% Remove a constraint
        function hobj = rmconstraint(hobj, ncon)
            % Remove one or several constraint strings
            if isnumeric(ncon)
                hobj.constraints = hobj.constraints(setdiff(1:end,ncon));
            end
        end
        
        %% Plot fit line
        function hobj = plot(hobj, varargin)
            % Plot result function
            if ishandle(hobj.fitline.handle), delete(hobj.fitline.handle); end
            hobj.fitline.handle = [];
            if ishandle(hobj.graphhandle)
                hobj.fitline.x = linspace(min(hobj.xdata),max(hobj.xdata),1000);                % x-values;
                hobj.fitline.y = hobj.fitfunction.call(hobj.parameters.values, hobj.fitline.x); % y-values
                hold(hobj.graphhandle,'on');
                % plot line
                hobj.fitline.handle = plot(hobj.graphhandle, hobj.fitline.x, hobj.fitline.y, 'tag', 'fitline');
                % set attributes
                for fn=fieldnames(hobj.fitline.plotstyle)'
                    set(hobj.fitline.handle, fn{1}, hobj.fitline.plotstyle.(fn{1})); 
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
                    bg.x = hobj.fitline.x;
                    bg.y = bgfunc(hobj.parameters.values(ismember(hobj.parameters.belongtocomponent,find(bgcomp))),hobj.fitline.x);
                    bg.handle = plot(hobj.graphhandle, bg.x, bg.y, 'tag', 'background');
                    bg.plotstyle = struct('color','k','linewidth',1,'linestyle','--');
                    hobj.graphelements{end+1} = bg;
                end
                    % ** Make a plot of more individual components ... (according "tag")
                % set attributes for all elements
                for iel=1:length(hobj.graphelements), 
                    for fn=fieldnames(hobj.graphelements{iel}.plotstyle)'
                        set(hobj.graphelements{iel}.handle, fn{1}, hobj.graphelements{iel}.plotstyle.(fn{1})); 
                    end
                end
            end    
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
