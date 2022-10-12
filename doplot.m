function doplot(plotstruct, opt)

% Low-level plotting function
% Content of plotstruct is plotted
% opt='keepview' holds range constant
%
% P. Steffens, 08/2014

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

%% merge datalist, if several slices (cell array)
% alldatalist is the single combined list for all plotting purposes in the following
% if iscell(plotstruct.datalist) 
%     if length(plotstruct.datalist)==1, alldatalist = plotstruct.datalist{1};
%     else alldatalist = cmbavg(plotstruct.datalist,'noAvg','ignore KF');
%     end
% else alldatalist = plotstruct.datalist;
% end
% if not, make cell array 
if ~iscell(plotstruct.datalist), m=plotstruct.datalist; plotstruct=rmfield(plotstruct,'datalist'); plotstruct.datalist{1}=m; clear m; end 
[ndata, ~, nbverts, pointinds, vertinds] = slicecounting(plotstruct.datalist);


%% Set linear or logarithmic scale

allcvalue = mergelist(plotstruct.datalist,'valuelist'); % merge datalist, if several slices (cell array)

if strcmpi(plotstruct.linlog,'LIN')
    cvalue = allcvalue(:,1);
elseif strcmpi(plotstruct.linlog,'LOG')
    m = min(allcvalue(:,1));
    cvalue = log10(allcvalue(:,1) + max(0,-m) + 1);
else
    fprintf('Error (in doplot): "lin" or "log" must be defined.\n');
    return;
end

%% Draw the coloured patch - loop over each slice

for nsl = 1:ndata
    
    switch upper(plotstruct.presentation)
        case 'VORONOI'
            if plotstruct.showcells, linestyle = '-'; else linestyle = 'none'; end
            patch('Faces',plotstruct.datalist{nsl}.faces + sum(nbverts(1:nsl-1)), ... 
                  'Vertices',plotstruct.vertexlist,'FaceColor','flat', ...
                  'FaceVertexCData',cvalue(pointinds{nsl}),'LineStyle',linestyle);
        case 'LINEAR'
     %       vertexcolor = griddatan(plotstruct.coordlist, cvalue, plotstruct.vertexlist, 'linear');
     %       patch('Faces',plotstruct.faces,'Vertices',plotstruct.vertexlist,'FaceColor','interp','FaceVertexCData',vertexcolor,'LineStyle','none');

            % ## Interpolation. ##
            % There are different ways to create an interpolatoin (see also options.m) in what concerns 
            % the algorithm, number of ppixels (if grid is used), etc. Apparently, the results can differ
            % significantly, and the parameters have to be set wisely.
            % Most rely on a Delaunay triangulation. Working with this has the advantage that one can put a 
            % size limit on triangles, which allows to limit the region of interpolation to areas that are not
            % too far away from the "nearest" real data point.
            % However, several possibilities:
            % This may be performed in the coordinate system of the dataset, or in the one of the plot. 
            % Both cases appear to make sense, but can give quite different results.
            % Implement both of them.

            % ## interpolation cannot be calculated for multiple slices ##
            % ** could rather easily do this on the individual ones ...


            [stdcell,stdratio] = getoption('stdcell','stdratio');

            switch upper(plotstruct.interpolationsystem)
                case 'DATA'  % work on plotstruct.datalist coordinates. Triangulation may already exist from elsewhere.

                    if isfield(plotstruct.datalist{nsl},'delaunaytri')
                        tri = plotstruct.datalist{nsl}.delaunaytri;          
                    else
                        stdratio = stdratio.(upper(plotstruct.datalist{nsl}.coordtype));
                        scalecoords = repmat(stdratio(:)', size(plotstruct.datalist{nsl}.coordlist,1), 1) .* plotstruct.datalist{nsl}.coordlist;
                        tri = delaunayn(scalecoords,{'QJ','Qbb','Qc'});
                    end
                    maxsize = stdcell.(upper(plotstruct.datalist{nsl}.coordtype));
                    coord = plotstruct.datalist{nsl}.coordlist;

                case 'PLOT'  % work on plot coordinates. Create triangulation for this purpose. (May nevertheless exist from previous interpolation actions here.)

                    if isfield(plotstruct, 'interpolatedimage') && iscell(plotstruct, 'interpolatedimage') ...
                            && length(plotstruct.interpolatedimage)>=nsl && isfield(plotstruct.interpolatedimage{nsl},'plottri')
                        tri = plotstruct.interpolatedimage{nsl}.plottri;          
                        % if sth is changed, 'interpolatedimage' has to be deleted (in fcplot etc); so long, this is safe
                    else
                        stdratio = stdratio.(upper(plotstruct.type));
                        scalecoords = repmat(stdratio(:)', size(plotstruct.coordlist,1), 1) .* plotstruct.coordlist(pointinds{nsl},:);
                        tri = delaunayn(scalecoords,{'QJ','Qbb','Qc'});
                    end
                    maxsize = stdcell.(upper(plotstruct.type));
                    coord = plotstruct.coordlist(pointinds{nsl},:);

                otherwise
                    fprintf('Interpolation options not recognized in plotstruct in doplot.\n');
            end

            % Now, mark the triangles that are too large. In those regions, no interpolation will be taken.
            for i=1:size(tri,1)  
                if    (max(coord(tri(i,:),1)) - min(coord(tri(i,:),1)))^2 + ...
                      (max(coord(tri(i,:),2)) - min(coord(tri(i,:),2)))^2 > plotstruct.interpolationlimit * sum(maxsize.^2)  
                                % this puts a limit on the size of the triangles (in relation to the maximum pixel size). 
                    toolarge(i) = true; %#ok<AGROW>
                else
                    toolarge(i) = false; %#ok<AGROW>
                end
            end
            toolarge = toolarge(1:i);
            % Plot the patch with interpolating colors between vertices
            switch upper(plotstruct.interpolationtype)
                case 'PATCHINTERP'
                    patch('Faces',tri(~toolarge,:),'Vertices',plotstruct.coordlist(pointinds{nsl},:), ...
                          'FaceColor','interp','FaceVertexCData',cvalue(pointinds{nsl}),'LineStyle','none');  % This is the interpolation provided by "patch"
                case {'PCOLORFACET','PCOLORSMOOTH'}  % Create our own interpolation, do everything explicitly by ourselves
                    ngr = plotstruct.interpolationgrid; % grid points in every dim.
                    if ~isfield(plotstruct,'interpolatedimage') % If not already calculated...
    %                     fprintf('Grid interpolation calculated based on the paramters from options.m\n');
                        [plotstruct.interpolatedimage{nsl}.xint,plotstruct.interpolatedimage{nsl}.yint] = ...
                            meshgrid(linspace(min(plotstruct.coordlist(pointinds{nsl},1)),max(plotstruct.coordlist(pointinds{nsl},1)),ngr), ...
                                     linspace(min(plotstruct.coordlist(pointinds{nsl},2)),max(plotstruct.coordlist(pointinds{nsl},2)),ngr));
                            % look wich points are enclosed in a triangle (too large ones have been thrown out above)
                        triind = tsearchn(plotstruct.coordlist(pointinds{nsl},1:2),tri,[plotstruct.interpolatedimage{nsl}.xint(:),plotstruct.interpolatedimage{nsl}.yint(:)]); 
                        plotstruct.interpolatedimage{nsl}.cvalint = nan(numel(plotstruct.interpolatedimage{nsl}.xint(:)),1);
                        triind(ismember(triind,find(toolarge))) = nan; % if point in a too large triangle, set nan
                        plotstruct.interpolatedimage{nsl}.cvalint(~isnan(triind)) = ... 
                            griddata(plotstruct.coordlist(pointinds{nsl},1) ,  plotstruct.coordlist(pointinds{nsl},2), ...
                                     cvalue(pointinds{nsl}), ...
                                     plotstruct.interpolatedimage{nsl}.xint(~isnan(triind)), plotstruct.interpolatedimage{nsl}.yint(~isnan(triind)), ...
                                     plotstruct.interpolationalgorithm);
        %                 connect=zeros((ngr-1)^2,4); % Faces (rectangles defined by four vertices)
        %                 for i=0:(ngr-2), connect(i*(ngr-1)+(1:(ngr-1)), 1:4) = [(1:(ngr-1))', (2:ngr)', ngr+(2:ngr)', ngr+(1:(ngr-1))']+ i*ngr; end
        %                 patch('Faces',connect,'Vertices',[xint(:),yint(:)],'FaceColor','interp','FaceVertexCData',cvalint,'LineStyle','none');
        %                 % 3 preceding lines appear equivalent to pcolor
                        plotstruct.interpolatedimage{nsl}.cvalint = reshape(plotstruct.interpolatedimage{nsl}.cvalint,size(plotstruct.interpolatedimage{nsl}.xint));
                        if strcmpi(plotstruct.interpolationsystem,'plot'), plotstruct.interpolatedimage{nsl}.plottri = tri; end % save tri
                    end
                    if strcmpi(plotstruct.interpolationtype,'pcolorsmooth')
                        phandle = pcolor(plotstruct.interpolatedimage{nsl}.xint,plotstruct.interpolatedimage{nsl}.yint,plotstruct.interpolatedimage{nsl}.cvalint);
                        shading interp; 
                    else
                        % little workaround for pcolor problem (not using full C-matrix on flat pixels)
                        xintstep = plotstruct.interpolatedimage{nsl}.xint(1,2)-plotstruct.interpolatedimage{nsl}.xint(1,1); 
                        yintstep = plotstruct.interpolatedimage{nsl}.yint(2,1)-plotstruct.interpolatedimage{nsl}.yint(1,1);                     
                        % Shift xint, yint half a pixel and extend matrices
                        phandle = pcolor([plotstruct.interpolatedimage{nsl}.xint([1,1:end],:)-xintstep/2, plotstruct.interpolatedimage{nsl}.xint([1,1:end],end)+xintstep/2], ...
                                         [plotstruct.interpolatedimage{nsl}.yint(:,[1,1:end])-yintstep/2; plotstruct.interpolatedimage{nsl}.yint(end,[1,1:end])+yintstep/2], ...
                                         [plotstruct.interpolatedimage{nsl}.cvalint,nan(ngr,1);nan(1,ngr+1)]);
                        shading flat
                    end
                    set(phandle,'linestyle','none')

                otherwise
                    fprintf('Interpolation options not recognized in plotstruct in doplot.\n');
            end
            % Show cell boundaries, if desired 
            if plotstruct.showcells
                patch('Faces',plotstruct.datalist{nsl}.faces + sum(nbverts(1:nsl-1)), ...
                      'Vertices',plotstruct.vertexlist,'FaceColor','none','FaceVertexCData',cvalue(pointinds{nsl}),'LineStyle','-');
            end
        case 'CONTOURF'
            [x,y]=meshgrid( linspace(min(plotstruct.coordlist(pointinds{nsl},1)),max(plotstruct.coordlist(pointinds{nsl},1)),100), ...
                            linspace(min(plotstruct.coordlist(pointinds{nsl},2)),max(plotstruct.coordlist(pointinds{nsl},2)),100));
            color = griddata( plotstruct.coordlist(pointinds{nsl},1), plotstruct.coordlist(pointinds{nsl},2), cvalue(pointinds{nsl}), x, y, 'linear');
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
    
    % If desired, plot slice edges
    if plotstruct.showedges
        if ~isfield(plotstruct,'sliceedges') || length(plotstruct.sliceedges)<nsl || isempty(plotstruct.sliceedges{nsl})
            % Calculate the edge line
            plotstruct.sliceedges{nsl} = findcontour(plotstruct.datalist{nsl}.faces);
        end
        for ncont = 1:length(plotstruct.sliceedges{nsl})
            if size(plotstruct.vertexlist,2) == 3 
                line('xdata',plotstruct.vertexlist(vertinds{nsl}(plotstruct.sliceedges{nsl}{ncont}),1), ...
                     'ydata',plotstruct.vertexlist(vertinds{nsl}(plotstruct.sliceedges{nsl}{ncont}),2), ...
                     'zdata',plotstruct.vertexlist(vertinds{nsl}(plotstruct.sliceedges{nsl}{ncont}),3), 'color', 'k', 'linewidth', 2);
            else
                line('xdata',plotstruct.vertexlist(vertinds{nsl}(plotstruct.sliceedges{nsl}{ncont}),1), ...
                     'ydata',plotstruct.vertexlist(vertinds{nsl}(plotstruct.sliceedges{nsl}{ncont}),2), 'color', 'k', 'linewidth', 2);
            end
        end
    end
    
end %loop over slices

%% Axes labels and title

switch upper(plotstruct.type)
    case 'QXY'
        if ~keepview
            axis equal; 
            axis([floor(min(plotstruct.vertexlist(:,1))*2)/2,ceil(max(plotstruct.vertexlist(:,1))*2)/2,  ...
              floor(min(plotstruct.vertexlist(:,2))*2)/2,ceil(max(plotstruct.vertexlist(:,2))*2)/2]);
            
        end
        xlabel(['Q_x (' char(197) '^{-1})']);
        ylabel(['Q_y (' char(197) '^{-1})']);
        if plotstruct.showvectors && isfield(plotstruct,'sampleinfo') 
            decorateQplot(plotstruct.sampleinfo.ax, plotstruct.sampleinfo.bx, plotstruct.sampleinfo.lattice, plotstruct.qvert, 'arrows', plotstruct.axeshandle, keepview); end  
        if plotstruct.showgrid    && isfield(plotstruct,'sampleinfo') 
            decorateQplot(plotstruct.sampleinfo.ax, plotstruct.sampleinfo.bx, plotstruct.sampleinfo.lattice, plotstruct.qvert, 'grid',   plotstruct.axeshandle, keepview); end  
    case 'HKLVECTORS'
        % draw orienting vectors
        if isfield(plotstruct,'sampleinfo')
            UB = UBmatrix( plotstruct.sampleinfo.lattice, plotstruct.sampleinfo.ax, plotstruct.sampleinfo.bx );
            basis1q = UB * plotstruct.basisvectors.vector1(:); basis2q = UB * plotstruct.basisvectors.vector2(:);
           
            if plotstruct.showvectors && size(plotstruct.vertexlist,2)==2   % in 2d
                 end1 = [plotstruct.sampleinfo.ax * UB' * basis1q/(basis1q'*basis1q), plotstruct.sampleinfo.ax * UB' * basis2q/(basis2q'*basis2q)];
                end2 = [plotstruct.sampleinfo.bx * UB' * basis1q/(basis1q'*basis1q), plotstruct.sampleinfo.bx * UB' * basis2q/(basis2q'*basis2q)];
                quiver([0,0], [0,0], [end1(1),end2(1)], [end1(2),end2(2)], 0, 'k');
                text(end2(1)+.1*end1(1),end2(2), num2str(plotstruct.sampleinfo.bx,'[%g,%g,%g]'),'Fontname','Verdana'); 
                if sign(end1(1))==-1; al='right'; else al='left'; end
                text(1.1*end1(1),end1(2),  num2str(plotstruct.sampleinfo.ax,'[%g,%g,%g]'), 'Fontname','Verdana','horizontalalignment',al);
            end
            if ~keepview, axis auto; asp=daspect; nasp=sqrt([basis2q'*basis2q,basis1q'*basis1q]); daspect([nasp(1:2)/mean(nasp(1:2))*mean(asp(1:2)), asp(3)] ); end
        end
    case 'ANGLES'
        xlabel('Scattering Angle (in-plane)');
        ylabel('Sample rotation angle');
        if ~keepview, asp=daspect; nasp=[1,1]; daspect([nasp(1:2)/mean(nasp(1:2))*mean(asp(1:2)), asp(3)] ); axis auto; end
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

% axis manual; % freeze limits

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


