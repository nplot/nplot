function doplot(plotstruct, opt)

% Low-level plotting function
% Content of plotstruct is plotted
% opt='keepview' holds range constant
%
% P. Steffens, 04/2009

%% To retain some old values
linecallback = get(plotstruct.scanlinehandle,'ButtonDownFcn');
prlinecallback = get(plotstruct.projectlinehandle,'ButtonDownFcn');
pp1callback = get(plotstruct.pp1handle,'ButtonDownFcn');
pp2callback = get(plotstruct.pp2handle,'ButtonDownFcn');
if nargin>1 && strcmpi(opt,'keepview')
    xrange = get(plotstruct.axeshandle,'xlim');
    yrange = get(plotstruct.axeshandle,'ylim');
    keepview = true;
else keepview = false;
end

%% Set current axes
axes(plotstruct.axeshandle);
cla(plotstruct.axeshandle);
if keepview
    set(gca,'xlim',xrange);
    set(gca,'ylim',yrange);
end
hold on
box on

%% Set linear or logarithmic scale
if strcmpi(plotstruct.linlog,'LIN')
    cvalue = plotstruct.datalist.valuelist(:,1);
elseif strcmpi(plotstruct.linlog,'LOG')
    m = min(plotstruct.datalist.valuelist(:,1));
    cvalue = log10(plotstruct.datalist.valuelist(:,1) + max(0,-m) + 1);
else
    fprintf('Error (in doplot): "lin" or "log" must be defined.\n');
    return;
end

%% Draw the coloured patch
switch upper(plotstruct.presentation)
    case 'VORONOI'
        if plotstruct.showcells, linestyle = '-'; else linestyle = 'none'; end
        patch('Faces',plotstruct.datalist.faces,'Vertices',plotstruct.vertexlist,'FaceColor','flat','FaceVertexCData',cvalue,'LineStyle',linestyle);
    case 'LINEAR'
 %       vertexcolor = griddatan(plotstruct.coordlist, cvalue, plotstruct.vertexlist, 'linear');
 %       patch('Faces',plotstruct.faces,'Vertices',plotstruct.vertexlist,'FaceColor','interp','FaceVertexCData',vertexcolor,'LineStyle','none');
        if isfield(plotstruct.datalist,'delaunaytri')
            tri = plotstruct.datalist.delaunaytri;          % ** Here and ...
        else
            stdratio = getoption('stdratio'); stdratio = stdratio.(upper(plotstruct.datalist.coordtype));
            scalecoords = repmat(stdratio(:)', size(plotstruct.datalist.coordlist,1), 1) .* plotstruct.datalist.coordlist;
            tri = delaunayn(scalecoords,{'QJ','Qbb','Qc'});   % ... here one might want to choose the qx,qy-system.
                                                            % In that case, one would need maintain an additional triangul. for this. 
                                                            % (The interpolated picture looks quite different depending on the system.) 
                                                            % ** Check the sense of this...!
        end
        i=1;
        stdcell = getoption('stdcell');
        maxsize = stdcell.(upper(plotstruct.datalist.coordtype));
        while i<=size(tri,1)  % delete triangles that are too large
            if    (max(plotstruct.coordlist(tri(i,:),1)) - min(plotstruct.coordlist(tri(i,:),1)))^2 + ...
                  (max(plotstruct.coordlist(tri(i,:),2)) - min(plotstruct.coordlist(tri(i,:),2)))^2 > 2*sum(maxsize.^2)  % ** this is arbitrary...
                tri = tri([1:(i-1),(i+1):end],:);
            else
                i=i+1;
            end
        end
        % Plot the patch with interpolating colors between vertices
        patch('Faces',tri,'Vertices',plotstruct.coordlist,'FaceColor','interp','FaceVertexCData',cvalue,'LineStyle','none');
        % Show cell boundaries, if desired 
        if plotstruct.showcells
            patch('Faces',plotstruct.datalist.faces,'Vertices',plotstruct.vertexlist,'FaceColor','none','FaceVertexCData',cvalue,'LineStyle','-');
        end
    case 'CONTOURF'
        [x,y]=meshgrid(linspace(min(plotstruct.coordlist(:,1)),max(plotstruct.coordlist(:,1)),100), linspace(min(plotstruct.coordlist(:,2)),max(plotstruct.coordlist(:,2)),100));
        color = griddata( plotstruct.coordlist(:,1), plotstruct.coordlist(:,2), cvalue, x, y, 'linear');
        contourf(x, y, color, 30, 'Linestyle', 'none');
    otherwise
        fprintf('Plot representation type not recognized in plotstruct in doplot.\n');
end
% If desired, plot dots at the positions of the measurement coordinates
if plotstruct.showcoordpoints
    if size(plotstruct.vertexlist,2)==3
        plot3(plotstruct.coordlist(:,1), plotstruct.coordlist(:,2), plotstruct.coordlist(:,3), '.k'); 
    else
        plot(plotstruct.coordlist(:,1), plotstruct.coordlist(:,2), '.k'); 
    end
end

%% Axes labels and title

switch upper(plotstruct.type)
    case 'QXY'
        if ~keepview
            axis([floor(min(plotstruct.vertexlist(:,1))*2)/2,ceil(max(plotstruct.vertexlist(:,1))*2)/2,  ...
              floor(min(plotstruct.vertexlist(:,2))*2)/2,ceil(max(plotstruct.vertexlist(:,2))*2)/2]);
%             axis equal; 
        end
        xlabel(['Q_x (' char(197) '^{-1})']);
        ylabel(['Q_y (' char(197) '^{-1})']);
        if plotstruct.showvectors  decorateQplot(plotstruct.datalist.sampleinfo.ax, plotstruct.datalist.sampleinfo.bx, plotstruct.datalist.sampleinfo.lattice, plotstruct.qvert, 'arrows', plotstruct.axeshandle, keepview); end  %#ok<SEPEX>
        if plotstruct.showgrid     decorateQplot(plotstruct.datalist.sampleinfo.ax, plotstruct.datalist.sampleinfo.bx, plotstruct.datalist.sampleinfo.lattice, plotstruct.qvert, 'grid',   plotstruct.axeshandle, keepview); end  %#ok<SEPEX>
    case 'HKLVECTORS'
        UB = UBmatrix( plotstruct.datalist.sampleinfo.lattice, plotstruct.datalist.sampleinfo.ax, plotstruct.datalist.sampleinfo.bx );
        basis1q = UB * plotstruct.basisvectors.vector1(:); basis2q = UB * plotstruct.basisvectors.vector2(:);
        % draw orienting vectors
        if plotstruct.showvectors && size(plotstruct.vertexlist,2)==2  % in 2d
            end1 = [plotstruct.datalist.sampleinfo.ax * UB' * basis1q/(basis1q'*basis1q), plotstruct.datalist.sampleinfo.ax * UB' * basis2q/(basis2q'*basis2q)];
            end2 = [plotstruct.datalist.sampleinfo.bx * UB' * basis1q/(basis1q'*basis1q), plotstruct.datalist.sampleinfo.bx * UB' * basis2q/(basis2q'*basis2q)];
            quiver([0,0], [0,0], [end1(1),end2(1)], [end1(2),end2(2)], 0, 'k');
            text(end2(1)+.1*end1(1),end2(2), num2str(plotstruct.datalist.sampleinfo.bx,'[%g,%g,%g]'),'Fontname','Verdana'); 
            if sign(end1(1))==-1; al='right'; else al='left'; end
            text(1.1*end1(1),end1(2),  num2str(plotstruct.datalist.sampleinfo.ax,'[%g,%g,%g]'), 'Fontname','Verdana','horizontalalignment',al);
        end
%         if ~keepview, asp=daspect; nasp=sqrt([basis2q'*basis2q,basis1q'*basis1q]); daspect([nasp(1:2)/mean(nasp(1:2))*mean(asp(1:2)), asp(3)] ); end
    case 'ANGLES'
        xlabel('Scattering Angle (in-plane)');
        ylabel('Sample rotation angle');
        if ~keepview, axis auto; end
    case 'QEPLANE'
        xlabel(['|Q| (' char(197) '^{-1})']);
        ylabel('E (meV)');
    case 'QXQYEN'
        xlabel(['Q_x (' char(197) '^{-1})']);
        ylabel(['Q_y (' char(197) '^{-1})']);
        zlabel('E (meV)');
        box on
        grid on
    case 'QXYZ'
        xlabel(['Q_x (' char(197) '^{-1})']);
        ylabel(['Q_y (' char(197) '^{-1})']);
        zlabel(['Q_z (' char(197) '^{-1})']);
        box on
        grid on
end

if isfield(plotstruct,'axesnames') && ~isempty(plotstruct.axesnames)
    xlabel(plotstruct.axesnames{1});
    ylabel(plotstruct.axesnames{2});
    if length(plotstruct.axesnames)>2, zlabel(plotstruct.axesnames{3}); end
end


    
%%

% Set callback function of all plot contents to the one of the axes
set( get(plotstruct.axeshandle,'Children'), 'ButtonDownFcn', get(plotstruct.axeshandle, 'ButtondownFcn'));

% If a scan line existed before, create a new one
if ~isempty(plotstruct.scanlinehandle)
    plotstruct.scanlinehandle = line( 'XData', plotstruct.scandef.xdat, 'YData', plotstruct.scandef.ydat, 'Marker','p','color','r', 'ButtonDownFcn', linecallback);     
end
% Similar for projection
if ~isempty(plotstruct.projectlinehandle)
    plotstruct.projectlinehandle = line( 'XData', plotstruct.scandef.xdat, 'YData', plotstruct.scandef.ydat, 'Marker','none','color','r', 'ButtonDownFcn', prlinecallback); 
    framevec = plotstruct.scandef.framevec; %(Assign these variables only to make notation shorter in following lines)
    xdat = plotstruct.scandef.xdat; ydat = plotstruct.scandef.ydat;
    plotstruct.projectframehandle = line( 'XData', [xdat(1)+framevec(1), xdat(2)+framevec(1), xdat(2)-framevec(1), xdat(1)-framevec(1), xdat(1)+framevec(1)], ...
                                          'YData', [ydat(1)+framevec(2), ydat(2)+framevec(2), ydat(2)-framevec(2), ydat(1)-framevec(2), ydat(1)+framevec(2)], 'Marker','none','color','r');
    plotstruct.pp1handle = line( 'XData', plotstruct.scandef.xdat(2),             'YData', plotstruct.scandef.ydat(2)            , 'Marker','o','MarkerFaceColor','r','color','r', 'ButtonDownFcn', pp1callback);
    plotstruct.pp2handle = line( 'XData', plotstruct.scandef.xdat(1)+framevec(1), 'YData', plotstruct.scandef.ydat(1)+framevec(2), 'Marker','o','MarkerFaceColor','r','color','r', 'ButtonDownFcn', pp2callback);
end

% Plot evtl. additional lines (powder etc.)
if isfield(plotstruct, 'linedef')
    if ~iscell(plotstruct.linedef), linedef{1}=plotstruct.linedef; else linedef = plotstruct.linedef; end  % ensure cell array
    for nl = 1:length(linedef)
        patch('Vertices',linedef{nl}.points, 'Faces',linedef{nl}.connection, linedef{nl}.properties{:});
    end
end

% Plot evtl. planes (volume cuts etc.)
if isfield(plotstruct, 'planedef')
    if ~iscell(plotstruct.planedef), planedef{1}=plotstruct.planedef; else planedef = plotstruct.planedef; end  % ensure cell array
    for npl = 1:length(planedef)
        xl=xlim; yl=ylim; zl=zlim;
        corners = [xl(1),yl(1),zl(1);xl(1),yl(1),zl(2);xl(1),yl(2),zl(1);xl(1),yl(2),zl(2);xl(2),yl(1),zl(1);xl(2),yl(1),zl(2);xl(2),yl(2),zl(1);xl(2),yl(2),zl(2)];
        [planevectors, origin] = getplaneparameter(planedef{npl}.normal(:)', planedef{npl}.c);
        % intersection
        [cutvertices, cutorder] = createmesh (corners, 1:8, planevectors', origin);
        % transform projected ccordinates back to full system
        cutvertices = cutvertices * planevectors' + repmat(origin(:)',size(cutvertices,1),1);
        patch('Vertices',cutvertices, 'Faces',cutorder, planedef{npl}.properties{:});
    end
end


% Save eventual changes to plotstruct
guidata(gcf, plotstruct); 


