function plotstruct = makeplotstruct(datalist,plottype,varargin)

% Setup a structure that contains the plotting information
%
% P. Steffens, 08/2011


[hbar, mass_n, meVJ] = getoption('hbar', 'mass_n', 'meVJ');

%% Create vertices and faces, if necessary
if ~isfield(datalist,'vertexlist') || ~isfield(datalist,'faces')
    % Calculate the graphical representation of 2D data
    % (i.e. calculate faces and vertices to plot a coloured patch)
    [stdcell,stdratio] = getoption('stdcell','stdratio','check',varargin);
    
    maxsize = stdcell.(upper(datalist.coordtype));
    [datalist.vertexlist,datalist.faces] = makevoronoi(datalist.coordlist(:,1:2), [maxsize(1), maxsize(2)],stdratio.(upper(datalist.coordtype))); 
     % ** this has to be changed
end

% %% In case of three-dimensional representation, obtain other dimensions for vertices
% % The (nontrivial) problem here is that the vertices, calculated during the process of setting up the patch, are obtained 
% % in two dimensions only. The third dimension is known for the data points (datalist.coordlist), but NOT for the vertices 
% % of the patch. These have to be inferred afterwards by interpolation, fitting, or whatever.
% %  (Note that the whole framework of cell patches etc. is intrinsically based on
% %   2 dimensions. The third dimension is used only for plotting.) 
% 
% if strcmp(upper(datalist.coordtype),'ANGLESENERGY')
%     %datalist.vertexlist(:,3) = griddata(datalist.coordlist(:,1), datalist.coordlist(:,2), datalist.coordlist(:,3), datalist.vertexlist(:,1), datalist.vertexlist(:,2));
%     %Linear fit to obtain missing coordinates.. 
%     %Do this perhaps more intelligently..! ** This does NOT work for multiple separate slices. **
%     [a1,a2,b]=fitplane(datalist.coordlist(:,1), datalist.coordlist(:,2), datalist.coordlist(:,3));
%     datalist.vertexlist(:,3) = a1 * datalist.vertexlist(:,1) + a2 * datalist.vertexlist(:,2) + b;
% end

%% Setup a structure that contains all plot information

plotstruct = getoption('plotopt'); %Standard plot settings from options-file
% plotstruct.faces = datalist.faces;
% plotstruct.valuelist = datalist.valuelist;
plotstruct.datalist = datalist;
plotstruct.type = plottype;

if ~isfield(datalist,'sampleinfo')
%     plotstruct.ax = datalist.sampleinfo.ax;
%     plotstruct.bx = datalist.sampleinfo.bx;
%     plotstruct.lattice = datalist.sampleinfo.lattice;
% else
    plotstruct.showvectors = false;
    plotstruct.showgrid = false;
end

if isfield(datalist,'sectionlist'), plotstruct.sectionlist= datalist.sectionlist; end

plotstruct.figurehandle = [];
plotstruct.axeshandle   = [];
plotstruct.scanlinehandle = [];
plotstruct.projectlinehandle = [];   % Central line for projection
plotstruct.pp1handle = [];           % End point of central line, used to rotate
plotstruct.pp2handle = [];           % Corner of frame, used to resize
plotstruct.projectframehandle = [];  % Frame that defines the integration area
plotstruct.scanaxeshandle = [];

plotstruct.linedef =  {};

%% Coordinate transformation into the desired coord. system for plotting

ctr = coordtransform(datalist,plottype);
if isempty(ctr), fprintf('No plot - conversion into desired coordinate axes not successful.\n'); plotstruct = []; return; end  
plotstruct.vertexlist = ctr.vertexlist;
plotstruct.coordlist = ctr.coordlist;



%% Depending on plot type, some additional information

plotstruct.scaninfo = [];
if any(strcmpi(plottype,{'QPLANE','QXY','Angles'}))
    plotstruct.qvert = datalist.QVERT;
    plotstruct.scaninfo = ['Q^\perp = ' num2str(datalist.QVERT)];
    if ~strcmp(num2str(datalist.QVERT),'0'); plotstruct.scaninfo = [plotstruct.scaninfo ' '  char(197) '^{-1}']; end
    if isfield(datalist,'KI') && isfield(datalist,'KF')
        plotstruct.scaninfo = [plotstruct.scaninfo ', Energy = ' num2str((datalist.KI^2-datalist.KF^2)*1E20*hbar^2/2/mass_n*meVJ) ' meV' ];
    end
elseif strcmpi(plottype,'QEPLANE') || strcmpi(plottype,'QXQYEN')
    plotstruct.qvert = datalist.QVERT;
    plotstruct.scaninfo = ['Q^\perp = ' num2str(datalist.QVERT)];
elseif strcmpi(plottype,'QXYZ')
    plotstruct.scaninfo = ['Energy = ' num2str((datalist.KI^2-datalist.KF^2)*1E20*hbar^2/2/mass_n*meVJ) ' meV' ];
end

plotstruct.scaninfo2 = { datalist.expname, ['Data: ' datalist.dataname]};
if any(strcmpi(datalist.constants, 'TEMP'))
    plotstruct.scaninfo2 = {plotstruct.scaninfo2{:},  ['Temperature: ' num2str(datalist.TEMP,4) ' K']};
end

if isfield(ctr,'variables')
    plotstruct.axesnames = ctr.variables;
end
    
    
    
