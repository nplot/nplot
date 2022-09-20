function avg = plotmultiple(scanlist,varargin)

% function avg = plotmultiple(scanlist,varargin)
%
% Create a plot that may consist of several scan files. The created data structure is returned. 
% scanlist  : Array of scan structures or as filenames as string (multi-filename format allowed).
%             Give cell array to plot several independent slices (for 3d plots)
% varargin  : Any pair of codeword/value. In particular
%             'coordinates': which coordinates to use. Allowed values are
%                           'angles' - work with a4 and psi(a3) (2d)
%                           'anglesEnergy' - like 'angles' plus Energy (3d) [a4,en,psi]
%                           'a4energy' - twotheta and Energy (2d)
%                           'anglesQz' - 'angles' plus Qz (3d)
%                           'qxy' - use qx and qy (in-plane momentum transfer) (2d)
%                           'qxyz' - like  'qxy' plus Qz (3d)
%                           'qxqyen' - like 'qxy' plus Energy (3d)
%                           'energyproj' - abs. value of Q, and energy (2d)
%                           'linearQ' - (IMPS) Q-coordinate along chosen direction, and energy (2d)
%                           'direct' - read variables from scan file. Give option 'variables' (cell array)
%             'plottype' : Which axes to use for plotting
%                           'qplane' - Qx,Qy coordinates (rec. Angst.)
%                           'qeplane' - Q,E axes
%                           'qxqyen' - Qx,Qy,En axes (3d)
%                           'qxqyqz' - Qx,Qy,Qz axes (3d)
% An automatic guess for appropriate coordinates and plottype is performed.
% Allowed codewords include most variables from the options.m file, which
% can thus be overwritten from the command line.
% Switches for gonio behaviour:
%             '2dmode' : for data using Mad's 2d-mode (gonios not zeroed)
%             'nogoniomode': completely ignore gonios - i.e. they're assumed flat or not installed 
%                           (note that's not equivalent to '2dmode' for instance if a3p installed)
%
% Examples:
% plotmultiple 0512[10,11,14:16]
% data = plotmultiple('0512[10,11,14:16]');
% plotmultiple({'012345','0123[46,48]'},'plottype qxqyen')
% plotmultiple 012345 normalizeto M1 normval 1000 vanacorr 2 det_eff [1,1,.99,1.01,1, ...]


% P. Steffens 09/2020


avg = [];
if ~iscell(scanlist) % ensure cell array
    if scanlist(1)=='{' && scanlist(end)=='}' %#ok<ALIGN> % simplified cell array arg.
        scanlist = regexpmatch(scanlist(2:end-1),'\w*\[[\w\,\:]+\][\w\.]*|\w+[\w\.]*');
    else m = scanlist; clear scanlist; scanlist{1}=m; end % single slice
end 
numberofslices = length(scanlist);

%% Preliminary stuff, for setting or input parameters

listtype = readinput('coordinates',varargin);
plottype = readinput('plottype',varargin);


if isempty(listtype) && isempty(plottype)
    
    if length(scanlist)>1, fprintf('If using several slices, ''plottype'' must be specified.\n'); return; end
    
    % first, for compatibility with earlier version, check if arguments without codewords are present
    % 2D plots
    if any(strcmpi(varargin,'energyproj')), listtype = 'A4Energy'; plottype = 'qeplane'; end
    if any(strcmpi(varargin,'a4energy')),   listtype = 'A4Energy'; plottype = 'qeplane'; end
    if any(strcmpi(varargin,'qxy')),        listtype = 'Qxy';      plottype = 'Qxy'; end

    if any(strcmpi(varargin,'linearq')),    listtype = 'linearq';  plottype = 'qeplane'; end
    
    if any(strcmpi(varargin,'marmot')),     listtype = 'CHEI';     plottype = 'CHEI'; end

    % 3D plots
    if any(strcmpi(varargin,'energy3d')),     listtype = 'AnglesEnergy'; plottype = 'qxqyen'; end
    if any(strcmpi(varargin,'anglesenergy')), listtype = 'AnglesEnergy'; plottype = 'qxqyen'; end


    if any(strcmpi(varargin,'qplane')),     plottype = 'Qxy'; end
    if any(strcmpi(varargin,'qeplane')),    plottype = 'qeplane'; end
    if any(strcmpi(varargin,'qxqyen')),     plottype = 'qxqyen'; end
    if any(strcmpi(varargin,'qxqyqz')),     plottype = 'qxyz'; end
end


% "Translate"...
if any(strcmpi({'q3d','qxyz','qxqyqz'},plottype)), plottype = 'qxyz'; end
if any(strcmpi({'energy3d','qxye','qxqyen'},plottype)), plottype = 'qxqyen'; end
    
% if coordinates and plottype not given, try to determine what is appropriate
if isempty(listtype) &&  ~isempty(plottype)
    % if plottype given, infer
    if any(strcmpi({'qplane','qxy'},plottype)), listtype = 'Qxy'; end
    if any(strcmpi({'qeplane','qe','energyproj','a4energy'},plottype)), listtype = 'a4energy'; end
    if any(strcmpi({'energy3d','qxye','qxqyen'},plottype)), listtype = 'AnglesEnergy'; end
    if any(strcmpi({'q3d','qxyz','qxqyqz'},plottype)), listtype = 'AnglesQZ'; end
end

if isempty(listtype)
    % otherwise try to guess
    filelist = multifilename(scanlist{1});
    s = tasread(filelist{1},'download');
    if isempty(s), fprintf('Error: No data in first file. Exit.\n'); avg=[]; return; end
    if (isfield(s.STEPS,'EN') && s.STEPS.EN ~= 0) || isfield(s.STEPS,'KI') || ...
       isfield(s.STEPS,'EI') || isfield(s.STEPS,'EF') || isfield(s.STEPS,'KF') 
        listtype = 'a4energy';
        plottype = 'qeplane';
        fprintf('Plot type set to Q-E-plane (using abs. value of Q). You can change by using ''coordinates'' or ''plottype'' option.\n');
    elseif isfield(s.STEPS,'A3') || isfield(s.STEPS,'A3P') || isfield(s.STEPS,'PSI') || isfield(s.STEPS,'QH') 
        listtype = 'Angles';
        plottype = 'Qxy';
    else    
        fprintf('Error: No plot, because type of plot could not be determined. Please use the options ''coordinates'' and ''plottype''.\nType ''help plotmultiple'' for help.\n');
        return;
    end
    
elseif isempty(plottype)  % listtype known, but not plottype
    if any(strcmpi(listtype,{'angles','qxy'})), plottype = 'Qxy'; end
    if any(strcmpi(listtype,{'energyproj','a4energy','linearq','qeplane'})), plottype = 'qeplane'; end
    if any(strcmpi(listtype,{'anglesEnergy','qxqyen'})), plottype = 'Energy3d'; end
    if any(strcmpi(listtype,{'AnglesQZ','qxyz'})), plottype = 'qxyz'; end
    
end
    

if isempty(listtype) || isempty(plottype)
    fprintf('Error: No plot, because type of plot could not be determined. Please use the options ''coordinates'' and ''plottype''.\nType ''help plotmultiple'' for help.\n');
    return;
end


% Some more command line options
if ~isempty(readinput('monitor',varargin)) && isempty(readinput('normalizeto',varargin)), val = readinput('monitor',varargin,'last');  varargin = [varargin, 'normalizeto M1','normval', val]; end
if ~isempty(readinput('time',varargin))    && isempty(readinput('normalizeto',varargin)), val = readinput('time',   varargin,'last');  varargin = [varargin, 'normalizeto TIME','normval', val]; end



%%
% Something to treat future case of several analyzers per channel 
% (e.g. Berlin-Multiflex, Panda-Bambus, camea, etc.). So far, only for Multiflexx
multianaoption = readinput('analyzers',varargin);
if ~isempty(multianaoption)
    if strcmpi(multianaoption,'all'), anachannels=1:5;
    elseif isnumeric(multianaoption), anachannels=multianaoption;
    else fprintf('Error: Could not read analyzers option. Exit.\n'); return; 
    end
    fprintf('Test-version of plotmultiple for multiflexx - use analyzers %s.\n',num2str(anachannels));
    numberofslices = numel(anachannels) * numberofslices; % each analyzer set becomes a slice
end


%% Starts here...

[stdcell,stdratio,stdbindist] = getoption('stdcell','stdratio','stdbindist','check',varargin);
% stdgrid = stdgrid.(upper(listtype));
stdbindist = stdbindist.(upper(listtype));
hassomedata = false;


% Loop over slices:
for number = 1:numberofslices

    % Takes scans (load, if necessary) and transform them into a linear list
    if isempty(multianaoption) %(this is the normal case)   
        linearlist = makelist(scanlist{number}, listtype, varargin{:});     
    else %(test version for several analyzer sets)
        linearlist = makelist(scanlist{mod(number-1,length(scanlist))+1}, listtype, varargin{:}, 'multiflexxanalyzer', anachannels(floor((number-1)/length(scanlist))+1));
    end

    % Combine and average data
    avg{number} = cmbavg(linearlist, [], varargin{:}); %#ok<*AGROW>
    if isempty(avg{number}); fprintf('No data in slice %d.\n',number);  continue; else hassomedata=true; end

    % For polarized data:
    if avg{number}.polarized
        calcstring = upper(readinput('calc',varargin));
        if isempty(calcstring) 
            if any(strcmpi('noplot',varargin)), continue; end
            fprintf('Error: For polarized data, select a pal-state or combination (use option ''calc'').\nCannot continue with multiple pal-states.\n');  if nargout < 1, clear avg; end, return;
        else [avg{number}, errorstate] = calclinearcombination(avg{number}, calcstring, 'output', true); if errorstate, if nargout < 1, clear avg; end, return; end
        end
    end
    
%     if any(strcmpi('noplot',varargin)), continue; end
    
    % Calculate the Voronoi cells
    ndims = size(avg{number}.coordlist,2);
    maxsize = stdcell.(upper(avg{number}.coordtype));   
    ratios  = stdratio.(upper(avg{number}.coordtype));
    constdims = (max(avg{number}.coordlist,[],1)- min(avg{number}.coordlist,[],1)) < stdbindist(:)';
            % Along these dimensions there is no variation. 
    if  ndims == 2 
        % Two-dimensional data, --> this is easy. The Voronoi cells are the final mesh to be plot
        [avg{number}.vertexlist, avg{number}.faces, avg{number}.delaunaytri] = ...
                makevoronoi(avg{number}.coordlist, maxsize, ratios, 'calctri');
        
    elseif ndims == 3 && any(constdims)
        % Three-dimensional, but one dimension constant, this is also easy, treat like 2D
        constdim = find(constdims,1,'first');   % This dim. is constant
        vardims = [1:(constdim-1),(constdim+1):ndims];  % These dims are treated as non-constant
        [avg{number}.vertexlist(:,vardims), avg{number}.faces, avg{number}.delaunaytri] = ...
                makevoronoi(avg{number}.coordlist(:,vardims), maxsize(vardims), ratios(vardims), 'calctri');
        avg{number}.vertexlist(:,constdim) = mean(avg{number}.coordlist(:,constdim));

        % ** Modify this case:
        % - allow for higher dimensions as long as only two vary
        % - fit a plane first and test if in this plane (useful e.g. for inclined planes...), as below
        
    else
        
        fprintf('This is the high-dimensional case.\n');
        
        % Do the full N-dim. Voronoi cells and then cut a slice to draw
        [vertexlist,cells] = makevoronoi(avg{number}.coordlist, maxsize, ratios);
        
        % Idea: test if points in a single plane. If yes, a single call to createmesh. 
        % If no, collect for each point its neighbors (delaunay), fit a
        % plane locally and obtain one face each (cutpoints of polyhypercut).
        % Or: if an equation f(r)=0 can be found for the surface, a call to createintersection might be much better. 
        % ** (needs to be done)
        
        planenormal = zeros(0,ndims);
        distvec = zeros(size(avg{number}.coordlist));
        
        % Fit a set of ndim-2 hyperplanes; their intersection defines a 2D-plane
        for npl = 1:(ndims-2)
            [planenormal(npl,1:ndims),planec(npl)] = fithyperplane(avg{number}.coordlist, planenormal);
            
            % Calculate distances to this plane
            dist = zeros(size(avg{number}.coordlist,1),1) - planec(npl);
            for di=1:ndims, dist = dist + avg{number}.coordlist(:,di) * planenormal(npl,di); end
            
            % Calculate distance vector for each point
            for di=1:ndims, distvec(:,di) = distvec(:,di) + dist*planenormal(npl,di); end
        end
        
        % Check if any distance vector is larger than allowed (for each coordinate)
        oneplane = true;
        for di=1:ndims
            if any(abs(distvec(:,di)) > stdbindist(di)), oneplane = false; end
        end
        
        if oneplane  % A single 2dim plane is to be used
            [planevectors, origin] = getplaneparameter(planenormal, planec); 
            planevectors = planevectors';
            [vertices, avg{number}.faces, numbers] = createmesh (vertexlist, cells, planevectors, origin);
            avg{number}.vertexlist = zeros(size(vertices,1),ndims);
            for i = 1:ndims % obtain coordinates in original system by linear combination of plane vectors
                avg{number}.vertexlist(:,i) = origin(i) + vertices(:,1)*planevectors(1,i) + vertices(:,2)*planevectors(2,i);
            end
            % Retain only those points that are relevant... (** do we want this? Or keep all for later?)
            avg{number}.coordlist = avg{number}.coordlist(numbers,:);
            avg{number}.valuelist = avg{number}.valuelist(numbers,:);
            avg{number}.monitorlist = avg{number}.monitorlist(numbers,:);
        else
            %...
            
            fprintf(['Abort: Was not able to find a hyperplane describing the data well enough.\nThis case is not yet implemented in plotmultiple.\n', ...
                     '(This may be a problem of too small tolerance; check stdbindist parameter in your options.m)\n']);
                 
            
        end
    end
    
%     %%**
%     avg{number}.sectionlist = number * ones(size(avg{number}.coordlist(:,1)));

end


if ~hassomedata
    fprintf('No data to plot.\n');  
else
    
%     avg = cmbavg(avg,'noAvg');  % **!! combine into single data structure?? No!
    if iscell(avg) && length(avg)==1, avg = avg{1}; end %** reduce if only one slice
    
    if isempty(avg), fprintf('No data to plot.\n'); 
    else
        % Plot
        if ~any(strcmpi('noplot',varargin)), fcplot (avg, plottype, varargin{:}); end
    end
end


if nargout < 1, clear avg; end



