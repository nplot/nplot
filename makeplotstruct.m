function plotstruct = makeplotstruct(datalist,plottype,varargin)

% Set up a structure that contains the plotting information

% P. Steffens, 08/2011 - 08/2014


[hbar, mass_n, meVJ] = getoption('hbar', 'mass_n', 'meVJ');

% if not, make cell array 
if ~iscell(datalist), m=datalist; clear datalist; datalist{1}=m; clear m; end 


%% Create vertices and faces, if necessary

if ~isfield(datalist{1},'vertexlist') || ~isfield(datalist{1},'faces')
    % Calculate the graphical representation of 2D data
    % (i.e. calculate faces and vertices to plot a coloured patch)
    [stdcell,stdratio] = getoption('stdcell','stdratio','check',varargin);
    
    maxsize = stdcell.(upper(datalist{1}.coordtype));
    [datalist{1}.vertexlist,datalist{1}.faces] = makevoronoi(datalist{1}.coordlist(:,1:2), [maxsize(1), maxsize(2)],stdratio.(upper(datalist{1}.coordtype))); 
     % ** this has to be changed (see plotmultiple)
     % ** it is here only for cases where fcplot is called from outside plotmultiple
end
% 
% % %% In case of three-dimensional representation, obtain other dimensions for vertices
% % % The (nontrivial) problem here is that the vertices, calculated during the process of setting up the patch, are obtained 
% % % in two dimensions only. The third dimension is known for the data points (datalist.coordlist), but NOT for the vertices 
% % % of the patch. These have to be inferred afterwards by interpolation, fitting, or whatever.
% % %  (Note that the whole framework of cell patches etc. is intrinsically based on
% % %   2 dimensions. The third dimension is used only for plotting.) 
% % 
% % if strcmp(upper(datalist.coordtype),'ANGLESENERGY')
% %     %datalist.vertexlist(:,3) = griddata(datalist.coordlist(:,1), datalist.coordlist(:,2), datalist.coordlist(:,3), datalist.vertexlist(:,1), datalist.vertexlist(:,2));
% %     %Linear fit to obtain missing coordinates.. 
% %     %Do this perhaps more intelligently..! ** This does NOT work for multiple separate slices. **
% %     [a1,a2,b]=fitplane(datalist.coordlist(:,1), datalist.coordlist(:,2), datalist.coordlist(:,3));
% %     datalist.vertexlist(:,3) = a1 * datalist.vertexlist(:,1) + a2 * datalist.vertexlist(:,2) + b;
% % end
% 
% 

%% Setup a structure that contains all plot information

plotstruct = getoption('plotopt','check',varargin); %Standard plot settings from options-file
% plotstruct.faces = datalist.faces;
% plotstruct.valuelist = datalist.valuelist;

plotstruct.datalist = datalist;     % always cell array now! changed from: [Use original "datalist" (may be cell array) (** ensure cell array !?!? **)]
plotstruct.type = plottype;



%% Coordinate transformation into the desired coord. system for plotting

for dnum=1:length(datalist)
    ctr{dnum} = coordtransform(datalist{dnum},plottype); %#ok<AGROW>
    if isempty(ctr), fprintf('No plot - conversion into desired coordinate axes not successful for slice %d.\n', dnum); plotstruct = []; return; end  
end
ctr = cmbavg(ctr,'noAvg','ignore KF');
plotstruct.vertexlist = ctr.vertexlist;
plotstruct.coordlist = ctr.coordlist;

%% initialize fields of plotstruct

if ~isfield(ctr,'sampleinfo')
    plotstruct.showvectors = false;
    plotstruct.showgrid = false;
else plotstruct.sampleinfo = ctr.sampleinfo;
end

if ~isfield(plotstruct,'interpolationtype')
    % only for compatibility with old options.m
    fprintf('Warning: missing entries in options.m file (probably using old version). Use defaults.\n');
    plotstruct.interpolationtype = 'pcolorsmooth';
    plotstruct.interpolationgrid = 200; 
    plotstruct.interpolationalgorithm = 'natural';
end

if isfield(ctr,'properties'), plotstruct.properties = ctr.properties; end

if isfield(ctr,'sectionlist'), plotstruct.sectionlist= ctr.sectionlist; end

plotstruct.figurehandle = [];
plotstruct.axeshandle   = [];
plotstruct.scanlinehandle = [];
plotstruct.projectlinehandle = [];   % Central line for projection
plotstruct.pp1handle = [];           % End point of central line, used to rotate
plotstruct.pp2handle = [];           % Corner of frame, used to resize
plotstruct.projectframehandle = [];  % Frame that defines the integration area
plotstruct.scanaxeshandle = [];

plotstruct.linedef =  {};

if iscell(plotstruct.datalist), plotstruct.sliceselection = true(1,length(plotstruct.datalist));
else plotstruct.sliceselection = true; end

%% Depending on plot type, some additional information

plotstruct.scaninfo = [];
if any(strcmpi(plottype,{'QPLANE','QXY','Angles'}))
    plotstruct.qvert = ctr.QVERT;
    plotstruct.scaninfo = ['Q^\perp = ' num2str(ctr.QVERT)];
    if ~strcmp(num2str(ctr.QVERT),'0'); plotstruct.scaninfo = [plotstruct.scaninfo ' '  char(197) '^{-1}']; end
    if isfield(ctr,'KI') && isfield(ctr,'KF')
        plotstruct.scaninfo = [plotstruct.scaninfo ', Energy = ' num2str((ctr.KI^2-ctr.KF^2)*1E20*hbar^2/2/mass_n*meVJ) ' meV' ];
    end
elseif strcmpi(plottype,'QEPLANE') || strcmpi(plottype,'QXQYEN')
    plotstruct.qvert = ctr.QVERT;
    plotstruct.scaninfo = ['Q^\perp = ' num2str(ctr.QVERT,'%4.2f')];
elseif strcmpi(plottype,'QXYZ')
    plotstruct.scaninfo = ['Energy = ' num2str((ctr.KI^2-ctr.KF^2)*1E20*hbar^2/2/mass_n*meVJ) ' meV' ];
end

plotstruct.scaninfo2 = { ctr.expname, ['Data: ' ctr.dataname]};
if any(strcmpi(ctr.constants, 'TEMP'))
    plotstruct.scaninfo2 = [plotstruct.scaninfo2,  ['Temperature: ' num2str(ctr.TEMP,'%5.1f') ' K']];
end

if isfield(ctr,'variables')
    plotstruct.axesnames = ctr.variables;
end
    
    
    
