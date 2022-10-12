function fh = fcplot(datalist,plottype,varargin)

% Color plot of the data in datalist
% (returns handle to the figure window)
%
% Plottype: 'qxy' (standard), 'angles', 'qeplane', 'qxqyen'

% P. Steffens 08/2011 - 09/2020

if nargin <2, plottype='QXY'; end

%% Create plot structure
plotstruct = makeplotstruct(datalist,plottype,varargin{:});
if isempty(plotstruct), return; end


%% Open new window and create contents
% fighandle = figure;
plotstruct.axeshandle = readinput('axeshandle',varargin);
newwindow = isempty(plotstruct.axeshandle);
if ~newwindow
    if ~ishandle(plotstruct.axeshandle) || ~strcmpi(get(plotstruct.axeshandle,'type'),'axes'), fprintf('Error: "axeshandle" is not valid.\n'); return; end
    plotstruct.figurehandle = get(plotstruct.axeshandle,'parent');
else
    plotstruct.figurehandle = figure('Position',[50 50 550 550],'Toolbar','figure');%,'Color',[1,1,1]);
    plotstruct.axeshandle   = axes('OuterPosition',[0 0 1 .95]);
end
fh = plotstruct.figurehandle;


% Test if 3D
if any(strcmpi(plotstruct.type,{'qxqyen','qxyz'})), threedim = true; else threedim = false; end


% Callback procedures, to be asociated to the menu and contextmenu entries
% (These are the "shorter" callbacks, they are simply defined as single
% strings. The longer callbacks are defined as separate procedures below.)
callback_linscale = ['plotstruct = guidata(gcf); plotstruct.linlog = ''LIN''; ' ...
    'if isfield(plotstruct,''interpolatedimage''), plotstruct=rmfield(plotstruct,''interpolatedimage''); end; guidata(gcf,plotstruct); doplot(plotstruct,''keepview''); clear plotstruct;' ]; 
callback_logscale = ['plotstruct = guidata(gcf); plotstruct.linlog = ''LOG''; ' ...
    'if isfield(plotstruct,''interpolatedimage''), plotstruct=rmfield(plotstruct,''interpolatedimage''); end; guidata(gcf,plotstruct); doplot(plotstruct,''keepview''); '... 
    'mrgl=mergelist(plotstruct.datalist,''valuelist''); if min(mrgl(:,1))<=0, fprintf(''Note: to plot on logarithmic scale, data are shifted to positive range by adding a constant.\n''); end; clear plotstruct; clear mrgl;']; 
callback_direct =   'plotstruct = guidata(gcf); plotstruct.presentation = ''voronoi''; guidata(gcf,plotstruct); doplot(plotstruct,''keepview''); clear plotstruct;' ; 
callback_interp =   'plotstruct = guidata(gcf); plotstruct.presentation = ''linear''; guidata(gcf,plotstruct); doplot(plotstruct,''keepview''); clear plotstruct;' ; 
callback_grid = 'plotstruct = guidata(gcf); plotstruct.showcells = ~plotstruct.showcells; guidata(gcf,plotstruct); doplot(plotstruct,''keepview''); clear plotstruct;' ;
callback_coord = 'plotstruct = guidata(gcf); plotstruct.showcoordpoints = ~plotstruct.showcoordpoints; guidata(gcf,plotstruct); doplot(plotstruct,''keepview''); clear plotstruct;'  ;
callback_vectors = 'plotstruct = guidata(gcf); plotstruct.showvectors = ~plotstruct.showvectors; guidata(gcf,plotstruct); doplot(plotstruct,''keepview''); clear plotstruct;' ;
callback_hklgrid = 'plotstruct = guidata(gcf); plotstruct.showgrid = ~plotstruct.showgrid; guidata(gcf,plotstruct); doplot(plotstruct,''keepview''); clear plotstruct;' ;
callback_edges = 'plotstruct = guidata(gcf); plotstruct.showedges = ~plotstruct.showedges; guidata(gcf,plotstruct); doplot(plotstruct,''keepview''); clear plotstruct;' ;
callback_save = ['plotstruct = guidata(gcf); [file,path] = uiputfile(''*.*'',''Save 2D Flatcone data''); ' ...
                 'if file~=0, sdat=[plotstruct.coordlist, mergelist(plotstruct.datalist,''valuelist'')]; save([path,file],''sdat'',''-ascii''); end; clear plotstruct;' ];

callback_scansave = ['plotstruct = guidata(gcf); [file,path] = uiputfile(''*.*'',''Save Flatcone scan data (simple format)''); ' ...
                     'if file~=0, sdat=[plotstruct.scandata.x, plotstruct.scandata.y, plotstruct.scandata.dy]; save([path,file],''sdat'',''-ascii''); end' ];
callback_np = ['plotstruct = guidata(gcf); np = inputdlg(''Enter number of scan points'',''Scan definition'',1,{num2str(plotstruct.scandef.np)}); ' ...
               'if ~isempty(np), plotstruct.scandef.np = str2num(np{1}); end; guidata(gcf, plotstruct); if isfield(plotstruct,''scandata''), doscanplot(plotstruct); end ' ];

callback_xax_auto =  'plotstruct = guidata(gcf); plotstruct.scandef.xaxiscoord = ''AUTO''; guidata(gcf,plotstruct); doscanplot(plotstruct);' ;
callback_xax_qm   =  'plotstruct = guidata(gcf); plotstruct.scandef.xaxiscoord = ''QMOD''; guidata(gcf,plotstruct); doscanplot(plotstruct);' ;
callback_xax_qx   =  'plotstruct = guidata(gcf); plotstruct.scandef.xaxiscoord = ''QX'';   guidata(gcf,plotstruct); doscanplot(plotstruct);' ;
callback_xax_qy   =  'plotstruct = guidata(gcf); plotstruct.scandef.xaxiscoord = ''QY'';   guidata(gcf,plotstruct); doscanplot(plotstruct);' ;
callback_xax_qh   =  'plotstruct = guidata(gcf); plotstruct.scandef.xaxiscoord = ''QH'';   guidata(gcf,plotstruct); doscanplot(plotstruct);' ;
callback_xax_qk   =  'plotstruct = guidata(gcf); plotstruct.scandef.xaxiscoord = ''QK'';   guidata(gcf,plotstruct); doscanplot(plotstruct);' ;
callback_xax_ql   =  'plotstruct = guidata(gcf); plotstruct.scandef.xaxiscoord = ''QL'';   guidata(gcf,plotstruct); doscanplot(plotstruct);' ;
callback_xax_qm   =  'plotstruct = guidata(gcf); plotstruct.scandef.xaxiscoord = ''QM'';   guidata(gcf,plotstruct); doscanplot(plotstruct);' ;
callback_xax_en   =  'plotstruct = guidata(gcf); plotstruct.scandef.xaxiscoord = ''EN'';   guidata(gcf,plotstruct); doscanplot(plotstruct);' ;

callcack_replotscan ='plotstruct = guidata(gcf); doscanplot(plotstruct);' ;

% Create menu and contextmenu (identical)
menuentry = uimenu('Label','FLATCONE');
contextmenu = uicontextmenu;
for f=[menuentry, contextmenu]
    uimenu(f,'Label','Logarithmic scale','Callback',callback_logscale);             % Change color scale to logarithmic
    uimenu(f,'Label','Linear scale','Callback',callback_linscale);                  % Change color scale to linear
    uimenu(f,'Label','Adjust Color Scale','Callback','colorlimits(gcbf);');         % Open window for coloraxis limits
    uimenu(f,'Label','Direct','Callback', callback_direct, 'Separator','on');       % Plot cells (voronoi)
    uimenu(f,'Label','Interpolated','Callback', callback_interp);                   % Interpolated plot
    uimenu(f,'Label','Set interpolation options','Callback',@setinterpolationoptions); % Dialog for interpolation options
    smen3 = uimenu(f,'Label','Graphical elements', 'Separator','on');
        uimenu(smen3,'Label','Data cells on/off','Callback', callback_grid);        %Show cell boundaries
        uimenu(smen3,'Label','Measured points on/off','Callback', callback_coord);  % Show coordinates of measured points
        uimenu(smen3,'Label','Orienting vectors on/off','Callback', callback_vectors);% Show vectors ax,bx
        uimenu(smen3,'Label','HKL-grid on/off','Callback', callback_hklgrid);       % Show underlying HKL-grid
        uimenu(smen3,'Label','Slice edges on/off','Callback', callback_edges);
    smen1 = uimenu(f,'Label','Caculate intersections');
        uimenu(smen1,'Label','Powder lines','Callback',@show_powder);               % Open dialog to choose powder lines to plot
        uimenu(smen1,'Label','Integer H,K,L planes','Callback',@callback_HKLintersect);
        uimenu(smen1,'Label','Delete all','Callback',@callback_deleteintersect);    % Remove all linedefs and planedefs
    uimenu(f,'Label','Delete data point by click', 'Callback', @start_clickdelete, 'Separator', 'on'); % Delete point by mouseclick
    uimenu(f,'Label','Select region to delete', 'Callback', @start_deleteregion);  % Select rectangular region to delete
    uimenu(f,'Label','Define a scan (interpolation)','Callback', @preparescan, 'Separator', 'on');  % Do the necessary to extract a 1D-Scan (interpolation)
    uimenu(f,'Label','Define a scan (integration)','Callback', @prepareintegration);                % Do the necessary to extract a 1D-Scan (integration)
    uimenu(f,'Label','Define a scan (projection)','Callback', @prepareprojection);                  % Do the necessary to extract a 1D-Scan (projection)
    smen4 = uimenu(f,'Label','Cut planes (in 3d)', 'Enable', 'off');
        uimenu(smen4,'Label','New (numerical input)','Callback',@addaplane);  
        uimenu(smen4,'Label','New vertical (x)','Callback',@addaplane);  
        uimenu(smen4,'Label','New vertical (y)','Callback',@addaplane);
        uimenu(smen4,'Label','Edit','Callback',[]);    
        uimenu(smen4,'Label','Delete','Callback',[]);  
        uimenu(smen4,'Label','Link to window (on/off)','Callback',[], 'Separator', 'on');  
        uimenu(smen4,'Label','Preview (on/off)','Callback',[]);  
        uimenu(smen4,'Label','Keep vertical (on/off)','Callback',[]);  
        uimenu(smen4,'Label','Do interpolation (per slice)','Callback',[], 'Separator', 'on');  
        uimenu(smen4,'Label','Do interpolation (global)','Callback',[]);
        uimenu(smen4,'Label','Do integration (per slice)','Callback',[]);
%       uimenu(smen4,'Label','Do integration (global)','Callback',[]);
        uimenu(smen4,'Label','Set parameters','Callback',[], 'Separator', 'on');
    uimenu(f,'Label','Add a vertical cut plane','Callback', @addaplane, 'Enable', 'off');           % Insert a cut plane
    uimenu(f,'Label','Save data to file', 'Callback', callback_save, 'Separator', 'on');        % Save 2D data to file
    smen2 = uimenu(f,'Label','Change axes', 'Separator', 'on');
        uimenu(smen2,'Label','Qx,Qy reciprocal Angstrom','Callback', @changeaxes);   % Change coordinate system
        uimenu(smen2,'Label','Angles 2Theta/Psi','Callback', @changeaxes);
        uimenu(smen2,'Label','Reciprocal lattice units','Callback', @changeaxes);
    
    if strcmpi(plotstruct.type,'qeplane')
        set(findobj('Label','Change axes', 'Parent', f), 'Enable', 'off');
        set(findobj('Label','HKL-grid on/off', 'Parent', smen3), 'Enable', 'off');
        set(findobj('Label','Orienting vectors on/off', 'Parent', smen3), 'Enable', 'off');
    end
    
    if threedim
%         set(findobj('Label','Interpolated', 'Parent', f), 'Enable', 'off');
%         set(findobj('Label','Delete data point by click', 'Parent', f), 'Enable', 'off');
%         set(findobj('Label','Select region to delete', 'Parent', f), 'Enable', 'off');
        set(findobj('Label','Orienting vectors on/off', 'Parent', smen3), 'Enable', 'off');
        set(findobj('Label','HKL-grid on/off', 'Parent', smen3), 'Enable', 'off');
%         set(findobj('Label','Define a scan (interpolation)', 'Parent', f), 'Enable', 'off');
%         set(findobj('Label','Define a scan (integration)', 'Parent', f), 'Enable', 'off');
        set(findobj('Label','Define a scan (projection)', 'Parent', f), 'Enable', 'off');
        set(findobj('Label','Add a vertical cut plane', 'Parent', f), 'Enable', 'on');
        set(findobj('Label','Change axes', 'Parent', f), 'Enable', 'off');
%         set(findobj('Label','Angles 2Theta/Psi', 'Parent', f), 'Enable', 'off');
    end
    
end
scsavemenu = [];
scexportmen = [];

set(plotstruct.axeshandle,'UIContextMenu',contextmenu);

% Callback to show mouse coordinates
if size(plotstruct.vertexlist,2)<3 % 2D-plot
    set(plotstruct.figurehandle, 'WindowButtonMotionFcn',  @showpos);
    plotstruct.poshandle = uicontrol('Style', 'Text', 'Units', 'Normalized', 'Position', [.65,0,.34,.05], 'String', '', ...
                                     'Horizontalalignment', 'right', 'BackgroundColor', get(plotstruct.figurehandle,'Color')); 
end



%%
% Save plotstruct as variable associated with this window
% (so that it can be retrieved later by guidata(..) for further use (replotting etc.)
guidata(plotstruct.figurehandle,plotstruct);


%% Plot
doplot(plotstruct);

% Title, Information line, etc.
title([' \fontsize{10}', plotstruct.scaninfo2, ['\fontsize{11} \bf ' plotstruct.scaninfo]]);
%uicontrol('Style','Text', 'Units', 'Normalized', 'Position',[0,.93,1,.07], ...
 %         'Tag', 'Titlebar2', 'String',plotstruct.scaninfo2);
% [normalizeto, normval] = getoption('normalizeto', 'normval', 'check', varargin);
if isfield(plotstruct,'properties') && isfield(plotstruct.properties,'normalization')
    annotation('Textbox','Units', 'pixels','Position',[1,1,550,25], 'String', ['Norm: ' plotstruct.properties.normalization], ...
           'Horizontalalignment', 'left', 'linestyle', 'none');
end

       
       
       
%% Callback to show mouse pointer coordinates

    function showpos(src,evnt)        
        if threedim, return; end
        try
        cp = get(plotstruct.axeshandle,'CurrentPoint');
        xrange = get(plotstruct.axeshandle,'xlim');
        yrange = get(plotstruct.axeshandle,'ylim');   
        
        if strcmpi(plotstruct.type, 'QEPLANE')
            if cp(1,1)>xrange(1) && cp(1,1)<xrange(2) && cp(1,2)>yrange(1) && cp(1,2)<yrange(2)
                set(plotstruct.poshandle, 'String', ['|q|  = ', num2str(cp(1,1),'%5.3f'), ',  E = ', num2str(cp(1,2),'%5.2f'), '    ']  );
            else 
                set(plotstruct.poshandle, 'String', '');
            end
            
        elseif strcmpi(plotstruct.type, 'QXY')
            if cp(1,1)>xrange(1) && cp(1,1)<xrange(2) && cp(1,2)>yrange(1) && cp(1,2)<yrange(2)
                UB = UBmatrix( plotstruct.sampleinfo.lattice, plotstruct.sampleinfo.ax, plotstruct.sampleinfo.bx);
                [H, K, L] = calcHKL( cp(1,1), cp(1,2), plotstruct.datalist{1}.QVERT, UB );
                set(plotstruct.poshandle, 'String', {['qxy = [', num2str(cp(1,1),'%5.3f'), ', ', num2str(cp(1,2),'%5.3f'), ']    '] , ...
                                                     ['HKL = ' num2str(H,'%6.3f') ', ' num2str(K,'%6.3f') ', ' num2str(L,'%6.3f') '    ' ] } );
            else 
                set(plotstruct.poshandle, 'String', '');
            end
            
        elseif strcmpi(plotstruct.type, 'HKLvectors')
            if cp(1,1)>xrange(1) && cp(1,1)<xrange(2) && cp(1,2)>yrange(1) && cp(1,2)<yrange(2)
                UB = UBmatrix( plotstruct.sampleinfo.lattice, plotstruct.sampleinfo.ax, plotstruct.sampleinfo.bx);
                hkl = plotstruct.basisvectors.origin + cp(1,1)*plotstruct.basisvectors.vector1 + cp(1,2)*plotstruct.basisvectors.vector2;
                qq = UB*hkl(:);
                set(plotstruct.poshandle, 'String', {['qxy = [', num2str(qq(1),'%5.3f'), ', ', num2str(qq(2),'%5.3f'), ']    '] , ...
                                                     ['HKL = ' num2str(hkl(1),'%6.3f') ', ' num2str(hkl(2),'%6.3f') ', ' num2str(hkl(3),'%6.3f') '    ' ] } );
            else 
                set(plotstruct.poshandle, 'String', '');
            end
            
        elseif strcmpi(plotstruct.type, 'ANGLES')
            if cp(1,1)>xrange(1) && cp(1,1)<xrange(2) && cp(1,2)>yrange(1) && cp(1,2)<yrange(2)
                set(plotstruct.poshandle, 'String', ['Angles = [', num2str(cp(1,1),'%5.2f'), ', ', num2str(cp(1,2),'%5.2f'), ']    '] );
            else 
                set(plotstruct.poshandle, 'String', '');
            end
        end
        catch
        end
    end

       
%% Callback routines for definition of scan line

   function mouseclick1(src,evnt)
     % Callback for Mouseclick in the 2D-Plot; starts defining a scan line
     if strcmp(get(gcf,'SelectionType'),'normal')
        plotstruct = guidata(gcf);
        set(gcf,'pointer','circle');
        cp = get(plotstruct.axeshandle,'CurrentPoint');
        xstart = cp(1,1);ystart = cp(1,2);
        plotstruct.scandef.xdat = [xstart,xstart];
        plotstruct.scandef.ydat = [ystart,ystart];
        if isempty(plotstruct.scanlinehandle), 
            plotstruct.scanlinehandle = line('XData',xstart,'YData',ystart,'Marker','p','color','r');
        else
            set(plotstruct.scanlinehandle,'XData',xstart,'YData',ystart,'Marker','p','color','r');
        end
        try %#ok<ALIGN>
        delete(plotstruct.projectlinehandle);   
        delete(plotstruct.pp1handle);          
        delete(plotstruct.pp2handle);          
        delete(plotstruct.projectframehandle); 
        catch end
        plotstruct.projectlinehandle=[]; plotstruct.pp1handle=[]; plotstruct.pp2handle=[]; plotstruct.projectframehandle=[];
        set(gcf,'WindowButtonMotionFcn',@mousemove1)
        set(gcf,'WindowButtonUpFcn',@mouseup1)
     end
 
        function mousemove1(src,evnt)
           showpos;
           cp = get(plotstruct.axeshandle,'CurrentPoint');
           plotstruct.scandef.xdat = [xstart,cp(1,1)];
           plotstruct.scandef.ydat = [ystart,cp(1,2)];
           set(plotstruct.scanlinehandle, 'XData', plotstruct.scandef.xdat, 'YData', plotstruct.scandef.ydat); drawnow expose
        end
   
        function mouseup1(src,evnt)
           % Reset behaviour of figure
           set(src,'Pointer','arrow')
           set(src,'WindowButtonMotionFcn',@showpos)
           set(src,'WindowButtonUpFcn','')
           % Set callback routine for line movement
           set(plotstruct.scanlinehandle, 'ButtonDownFcn', @mouseclick_line);
           %** {1}
           plotstruct.scandef.np = estimateNP([plotstruct.scandef.xdat', plotstruct.scandef.ydat'], plotstruct.vertexlist, plotstruct.datalist{1}.faces);
           guidata(gcf, plotstruct);
           doscanplot(plotstruct);
           set(scsavemenu, 'Enable', 'on');
           set(scexportmen, 'Enable', 'on');
        end
              
   end %mouseclick1

%% Callback routines for movement of scan line

    function mouseclick_line(src,event)
        if strcmp(get(gcf,'SelectionType'),'normal')
            plotstruct = guidata(gcf);
            % Get current mouse position
            cp = get(plotstruct.axeshandle,'CurrentPoint');
            xstart = cp(1,1);ystart = cp(1,2);
            oldpointer = get(gcf,'pointer');
            set(gcf,'pointer','fleur');
            set(gcf,'WindowButtonMotionFcn',@mousemove_sl)
            set(gcf,'WindowButtonUpFcn',@mouseup_sl)
        end
        
        function mousemove_sl(src,event)
           showpos;
           cp = get(plotstruct.axeshandle,'CurrentPoint');
           % Shift start and end points like mouse movement (relative to start values)
           plotstruct.scandef.xdat = plotstruct.scandef.xdat + cp(1,1) - xstart;
           plotstruct.scandef.ydat = plotstruct.scandef.ydat + cp(1,2) - ystart;
           xstart = cp(1,1); ystart = cp(1,2);
           set(plotstruct.scanlinehandle, 'XData', plotstruct.scandef.xdat, 'YData', plotstruct.scandef.ydat); drawnow
        end
        
        function mouseup_sl(src,evnt)
           set(src,'WindowButtonMotionFcn',@showpos)
           set(src,'WindowButtonUpFcn','')
           set(gcf,'pointer', oldpointer')
           guidata(gcf, plotstruct);
           doscanplot(plotstruct);
        end
        
    end %mouseclick_line


%% Callback routines for definition of projection

   function mouseclick2(src,evnt)
     % Callback for Mouseclick in the 2D-Plot; starts defining a projection line
     if strcmp(get(gcf,'SelectionType'),'normal')
        plotstruct = guidata(gcf);
        set(gcf,'pointer','circle');
        cp = get(plotstruct.axeshandle,'CurrentPoint');
        xstart = cp(1,1);ystart = cp(1,2);
        framevec = [cp(1,1), cp(1,2)];              
        try %#ok<ALIGN>
        delete(plotstruct.scanlinehandle);
        delete(plotstruct.projectlinehandle);   
        delete(plotstruct.pp1handle);          
        delete(plotstruct.pp2handle);          
        delete(plotstruct.projectframehandle); 
        catch end %#ok<*SEPEX>
        plotstruct.scanlinehandle=[];
        plotstruct.projectlinehandle = line('XData',xstart,'YData',ystart,'Marker','none','color','r');
        plotstruct.projectframehandle = line('XData',xstart,'YData',ystart,'Marker','none','color','r');
        plotstruct.pp1handle = line('XData',xstart,'YData',ystart,'Marker','o','MarkerFaceColor','r','color','r');
        plotstruct.pp2handle = line('XData',xstart,'YData',ystart,'Marker','o','MarkerFaceColor','r','color','r');
        set(gcf,'WindowButtonMotionFcn',@mousemove2)
        set(gcf,'WindowButtonUpFcn',@mouseup2)
     end
 
        function mousemove2(src,evnt)
           showpos;
           cp = get(plotstruct.axeshandle,'CurrentPoint');
           plotstruct.scandef.xdat = [xstart,cp(1,1)];
           plotstruct.scandef.ydat = [ystart,cp(1,2)];
           set(plotstruct.projectlinehandle, 'XData', plotstruct.scandef.xdat, 'YData', plotstruct.scandef.ydat); 
           yl = get(plotstruct.axeshandle,'ylim'); xl = get(plotstruct.axeshandle,'xlim');
           vec = [plotstruct.scandef.xdat(2)-plotstruct.scandef.xdat(1), plotstruct.scandef.ydat(2)-plotstruct.scandef.ydat(1)];
           framevec = [-vec(2) / (yl(2)-yl(1))*(xl(2)-xl(1)), vec(1) / (xl(2)-xl(1))*(yl(2)-yl(1))] * min(0.5 , 0.25*sqrt(sum([vec(1)/(xl(2)-xl(1)), vec(2)/(yl(2)-yl(1))].^2)) );
           plotstruct.scandef.framevec = framevec;
           set(plotstruct.projectframehandle, 'XData', [xstart+framevec(1), cp(1,1)+framevec(1), cp(1,1)-framevec(1), xstart-framevec(1), xstart+framevec(1)], ...
                                              'YData', [ystart+framevec(2), cp(1,2)+framevec(2), cp(1,2)-framevec(2), ystart-framevec(2), ystart+framevec(2)]);
           set(plotstruct.pp1handle, 'XData', cp(1,1), 'YData', cp(1,2));
           set(plotstruct.pp2handle, 'XData', xstart+framevec(1), 'YData', ystart+framevec(2));           
           drawnow
        end
   
        function mouseup2(src,evnt)
           % Reset behaviour of figure
           set(src,'Pointer','arrow')
           set(src,'WindowButtonMotionFcn',@showpos)
           set(src,'WindowButtonUpFcn','')
           % Set callback routine for line movement
           set(plotstruct.projectlinehandle, 'ButtonDownFcn', @mouseclick_projectline);
           set(plotstruct.pp1handle, 'ButtonDownFcn', @mouseclick_pp1);
           set(plotstruct.pp2handle, 'ButtonDownFcn', @mouseclick_pp2);
           % ** ({1})
           plotstruct.scandef.np = estimateNP([plotstruct.scandef.xdat', plotstruct.scandef.ydat'], plotstruct.vertexlist, plotstruct.datalist{1}.faces);
           guidata(gcf, plotstruct);
           doscanplot(plotstruct);
           set(scsavemenu, 'Enable', 'on');
           set(scexportmen, 'Enable', 'on');
        end
              
   end %mouseclick2

%% Callback routines for movement, rotation and resize of projection 

    function mouseclick_projectline(src,event) %#ok<*INUSD>
        if strcmp(get(gcf,'SelectionType'),'normal')
            plotstruct = guidata(gcf);
            % Get current mouse position
            cp = get(plotstruct.axeshandle,'CurrentPoint');
            xstart = cp(1,1);ystart = cp(1,2);
            oldpointer = get(gcf,'pointer');
            set(gcf,'pointer','fleur');
            set(gcf,'WindowButtonMotionFcn',@mousemove_pl)
            set(gcf,'WindowButtonUpFcn',@mouseup_pl)
        end
        
        function mousemove_pl(src,event)
           showpos;
           cp = get(plotstruct.axeshandle,'CurrentPoint');
           % Shift all points of projection definition according to mouse movement (relative to start values)
           plotstruct.scandef.xdat = plotstruct.scandef.xdat + cp(1,1) - xstart;
           plotstruct.scandef.ydat = plotstruct.scandef.ydat + cp(1,2) - ystart;
           xstart = cp(1,1); ystart = cp(1,2);
           set(plotstruct.projectlinehandle, 'XData', plotstruct.scandef.xdat, 'YData', plotstruct.scandef.ydat); 
           framevec = plotstruct.scandef.framevec; %(Assign these variables only to make notation shorter in following lines)
           xdat = plotstruct.scandef.xdat; ydat = plotstruct.scandef.ydat;
           set(plotstruct.projectframehandle, 'XData', [xdat(1)+framevec(1), xdat(2)+framevec(1), xdat(2)-framevec(1), xdat(1)-framevec(1), xdat(1)+framevec(1)], ...
                                              'YData', [ydat(1)+framevec(2), ydat(2)+framevec(2), ydat(2)-framevec(2), ydat(1)-framevec(2), ydat(1)+framevec(2)]);
           set(plotstruct.pp1handle, 'XData', plotstruct.scandef.xdat(2),             'YData', plotstruct.scandef.ydat(2));
           set(plotstruct.pp2handle, 'XData', plotstruct.scandef.xdat(1)+framevec(1), 'YData', plotstruct.scandef.ydat(1)+framevec(2));
           drawnow
        end
        
        function mouseup_pl(src,evnt)
           set(src,'WindowButtonMotionFcn',@showpos)
           set(src,'WindowButtonUpFcn','')
           set(gcf,'pointer',oldpointer)
           guidata(gcf, plotstruct);
           doscanplot(plotstruct);
        end
        
    end %mouseclick_projectline



    function mouseclick_pp1(src,event)  % Rotates the whole projection
        if strcmp(get(gcf,'SelectionType'),'normal')
            plotstruct = guidata(gcf);
            set(gcf,'WindowButtonMotionFcn',@mousemove_pp1)
            set(gcf,'WindowButtonUpFcn',@mouseup_pp1)
            % Get current mouse position
            cp = get(plotstruct.axeshandle,'CurrentPoint');
            xstart = cp(1,1);ystart = cp(1,2);
        end
        
        function mousemove_pp1(src,event)
           showpos;
           cp = get(plotstruct.axeshandle,'CurrentPoint');
           yl = get(plotstruct.axeshandle,'ylim'); xl = get(plotstruct.axeshandle,'xlim'); % to consider axis ratio
           % Shift all points of projection definition according to mouse movement (relative to start values)
           plotstruct.scandef.xdat(2) = cp(1,1);
           plotstruct.scandef.ydat(2) = cp(1,2);
           vec_alt = [xstart-plotstruct.scandef.xdat(1), ystart-plotstruct.scandef.ydat(1)] ./ [xl(2)-xl(1), yl(2)-yl(1)];
           xstart = cp(1,1); ystart = cp(1,2);
           vec_neu = [xstart-plotstruct.scandef.xdat(1), ystart-plotstruct.scandef.ydat(1)] ./ [xl(2)-xl(1), yl(2)-yl(1)];
           % Determine the angle of rotation of the central line
           cosa = vec_alt*vec_neu' / sqrt(sum(vec_alt.^2)) / sqrt(sum(vec_neu.^2));
           signa = sign(vec_neu(2)*vec_alt(1)-vec_neu(1)*vec_alt(2));
           % Rotate framevec accordingly
           framevec = plotstruct.scandef.framevec ./ [xl(2)-xl(1), yl(2)-yl(1)]  * [ cosa, signa*sqrt(1-cosa^2); -signa*sqrt(1-cosa^2), cosa] .* [xl(2)-xl(1), yl(2)-yl(1)];
           plotstruct.scandef.framevec = framevec;
           set(plotstruct.projectlinehandle, 'XData', plotstruct.scandef.xdat, 'YData', plotstruct.scandef.ydat); 
           xdat = plotstruct.scandef.xdat; ydat = plotstruct.scandef.ydat; %(Assign these variables only to make notation shorter in following lines)
           set(plotstruct.projectframehandle, 'XData', [xdat(1)+framevec(1), xdat(2)+framevec(1), xdat(2)-framevec(1), xdat(1)-framevec(1), xdat(1)+framevec(1)], ...
                                              'YData', [ydat(1)+framevec(2), ydat(2)+framevec(2), ydat(2)-framevec(2), ydat(1)-framevec(2), ydat(1)+framevec(2)]);
           set(plotstruct.pp1handle, 'XData', plotstruct.scandef.xdat(2),             'YData', plotstruct.scandef.ydat(2));
           set(plotstruct.pp2handle, 'XData', plotstruct.scandef.xdat(1)+framevec(1), 'YData', plotstruct.scandef.ydat(1)+framevec(2));
           drawnow
        end
        
        function mouseup_pp1(src,evnt)
           set(src,'WindowButtonMotionFcn',@showpos)
           set(src,'WindowButtonUpFcn','')
           guidata(gcf, plotstruct);
           doscanplot(plotstruct);
        end
        
    end %mouseclick_pp1


    function mouseclick_pp2(src,event)  % Changes the frame for the projection
        if strcmp(get(gcf,'SelectionType'),'normal')
            plotstruct = guidata(gcf);
            set(gcf,'WindowButtonMotionFcn',@mousemove_pp2)
            set(gcf,'WindowButtonUpFcn',@mouseup_pp2)
        end
        
        function mousemove_pp2(src,event)
           showpos;
           cp = get(plotstruct.axeshandle,'CurrentPoint');
           % Shift all points of projection definition according to mouse movement (relative to start values)
           framevec = [cp(1,1)-plotstruct.scandef.xdat(1), cp(1,2)-plotstruct.scandef.ydat(1)];
           plotstruct.scandef.framevec = framevec;
           set(plotstruct.projectlinehandle, 'XData', plotstruct.scandef.xdat, 'YData', plotstruct.scandef.ydat); 
           xdat = plotstruct.scandef.xdat; ydat = plotstruct.scandef.ydat; %(Assign these variables only to make notation shorter in following lines)
           set(plotstruct.projectframehandle, 'XData', [xdat(1)+framevec(1), xdat(2)+framevec(1), xdat(2)-framevec(1), xdat(1)-framevec(1), xdat(1)+framevec(1)], ...
                                              'YData', [ydat(1)+framevec(2), ydat(2)+framevec(2), ydat(2)-framevec(2), ydat(1)-framevec(2), ydat(1)+framevec(2)]);
           set(plotstruct.pp1handle, 'XData', plotstruct.scandef.xdat(2),             'YData', plotstruct.scandef.ydat(2));
           set(plotstruct.pp2handle, 'XData', plotstruct.scandef.xdat(1)+framevec(1), 'YData', plotstruct.scandef.ydat(1)+framevec(2));
           drawnow
        end
        
        function mouseup_pp2(src,evnt)
           set(src,'WindowButtonMotionFcn',@showpos)
           set(src,'WindowButtonUpFcn','')
           guidata(gcf, plotstruct);
           doscanplot(plotstruct);
        end
        
    end %mouseclick_pp1


%% Callback : Add a plane
    function addaplane(src,evnt)
        plotstruct = guidata(gcf);
        newplanedef.normal = [1;0;0];
        newplanedef.c = sum(xlim)/2;
        newplanedef.properties = {'FaceColor','r','FaceAlpha',.35,'Tag','cutplane','ButtonDownFcn',@mouseclickonplane}; %** use options.m for this?
        if ~isfield(plotstruct,'planedef'), plotstruct.planedef = {}; end
        plotstruct.planedef = [plotstruct.planedef, newplanedef];
        guidata(gcf, plotstruct);
        doplot(plotstruct);
    end %addaplane

%% Callback for plane movements

function mouseclickonplane(src,evnt)
     % Callback for Mouseclick on a plane
     if strcmp(get(gcf,'SelectionType'),'normal')
        plotstruct = guidata(gcf);
        cp = get(plotstruct.axeshandle,'CurrentPoint');
        camPos= get(plotstruct.axeshandle, 'CameraPosition');  
        xyzlim = axis; % axis limits
        corners = xyzlim([1,3,5;1,3,6;1,4,5;1,4,6;2,3,5;2,3,6;2,4,5;2,4,6]);
        axunits = get(gca,'units'); set(gca,'units','pixels'); axpos = get(gca,'position'); set(gca,'units',axunits);
        rot = rotationtoviewingplane(plotstruct.axeshandle);
        cornersr = rot*corners';
        pixelscale = [axpos(3)/(max(cornersr(1,:))-min(cornersr(1,:))), axpos(4)/(max(cornersr(2,:))-min(cornersr(2,:)))];
        pointh=[]; lineh=[]; %handles for graphics during movement
        
        % first, look on which plane is clicked
        mindist = inf; plane=nan;
        for npl = 1:length(plotstruct.planedef)
            % Point at which the vieweing line intersects
            pt = cp(1,:) + (plotstruct.planedef{npl}.c - cp(1,:)*plotstruct.planedef{npl}.normal)/((cp(2,:)-cp(1,:))*plotstruct.planedef{npl}.normal)*(cp(2,:)-cp(1,:));
            if pt(1)<xyzlim(1) || pt(1)>xyzlim(2) || pt(2)<xyzlim(3) || pt(2)>xyzlim(4) || pt(3)<xyzlim(5) || pt(3)>xyzlim(6)
                continue; end % if not in the visible range
            if norm(pt-camPos)<mindist, plane=npl; mindist=norm(pt-camPos); hitpoint=pt; end
        end
        if isnan(plane), return; end
        % Then, check if this place is "visible", i.e. no data in front of it
        hitface = findface3d (plotstruct.axeshandle, plotstruct.coordlist, plotstruct.vertexlist, plotstruct.datalist);
        if ~isnan(hitface) && norm(plotstruct.coordlist(hitface,:)-camPos)<mindist, return; end
        
        % get handle to this plane object (may be different from calling obj.)
        thisplaneh = findobj(get(plotstruct.axeshandle,'children'),'Tag','cutplane'); 
        thisplaneh = thisplaneh(end+1-plane);
        oldcolor = get(thisplaneh,'facecolor');
        set(thisplaneh,'facecolor','b');
        oldlimmode = get(plotstruct.axeshandle,'xlimmode');
        axis manual; % fix axes limits
        
        % Determine, if clicked on a border or somewhere in the middle
        % According to that, do perpendicular movement or rotation
        vertp = get(thisplaneh,'vertices');
        facp = get(thisplaneh,'faces');
        % determine distance of hitpoint to each border
        borderdist = nan(1,numel(facp));
        for nb= 1:numel(facp)
            gvec = vertp(facp(nb),:)-vertp(facp(mod(nb,numel(facp))+1),:); %connection of the two endpoints, describes border line
            distvec = (hitpoint-vertp(facp(nb),:)) - ((hitpoint-vertp(facp(nb),:))*gvec')/(gvec*gvec')*gvec;
            % this distance in pixels
            distvec = rot*distvec'; borderdist(nb) = sqrt((distvec(1)*pixelscale(1))^2+(distvec(2)*pixelscale(2))^2);
        end
        % determine distance of hitpoint to each corner
        cornerdist = nan(1,size(vertp,1));
        for nc= 1:size(vertp,1)
            % this distance in pixels
            distvec = rot*(hitpoint-vertp(nc,:))'; cornerdist(nc) = sqrt((distvec(1)*pixelscale(1))^2+(distvec(2)*pixelscale(2))^2);
        end
        
        if any(cornerdist<5)    % if clicked on a corner
            % move this cornerpoint, while keeping the two neighbors fixed (might also keep an edge fixed)
            [~,nc] = min(cornerdist);
            pointh(1) = plot3(vertp(nc,1),vertp(nc,2),vertp(nc,3),'.r','markersize',20);
            ncneighb = facp( [mod(find(facp==nc),numel(facp))+1, mod(find(facp==nc)-2,numel(facp))+1] ); % neighbors are right and left in facp
            pointh(2) = plot3(vertp(ncneighb(1),1),vertp(ncneighb(1),2),vertp(ncneighb(1),3),'.b','markersize',20);
            pointh(3) = plot3(vertp(ncneighb(2),1),vertp(ncneighb(2),2),vertp(ncneighb(2),3),'.b','markersize',20);
            % identify line on which to move (edge of the coord box)
            mvdim = find([~any(abs(vertp(nc,1)-xyzlim(1:2))<1e-9); ~any(abs(vertp(nc,2)-xyzlim(3:4))<1e-9); ~any(abs(vertp(nc,3)-xyzlim(5:6))<1e-9)]);            
            if isempty(mvdim), mvdim=1; end; mvdim=mvdim(1);
            set(gcf,'WindowButtonMotionFcn',@mousemovepl2);
            
        elseif any(borderdist<5)    % if clicked on a border
            % move this edge, while keeping one point fixed
            [~,nb] = min(borderdist); endp = [facp(nb), facp(mod(nb,numel(facp))+1)];
            lineh = line(vertp(endp,1),vertp(endp,2),vertp(endp,3),'color','r','linewidth',3);
            % point which is furthest away from this line
            dist = -inf(1,size(vertp,1));
            for nc = setdiff(1:size(vertp,1), endp)
                dist(nc) = norm(cross(vertp(nc,:)-vertp(endp(1),:),vertp(endp(2),:)-vertp(endp(1),:))); % (not true distance; scaled)
            end
            [~,nc] = max(dist);
            pointh = plot3(vertp(nc,1),vertp(nc,2),vertp(nc,3),'.b','markersize',20);
            % identify face on which the line lies
            cdim = find( abs(vertp(endp(1),:)-vertp(endp(2),:))<1e-9 & ...  %(const. dim)
                        [any(abs(vertp(endp(1),1)-xyzlim(1:2))<1e-9), any(abs(vertp(endp(1),2)==xyzlim(3:4))<1e-9), any(abs(vertp(endp(1),3)==xyzlim(5:6))<1e-9)] ,1); 
            axface2d = sortrows(unique(corners(:,setdiff(1:3,cdim)),'rows'),[1,2]); axface2d=axface2d([1,2,4,3],:); % right order!
            dvec=vertp(endp(2),:)-vertp(endp(1),:); %nvec2d=dvec(setdiff(1:3,cdim)); %nvec2d=[-nvec2d(2),nvec2d(1)];
            set(gcf,'WindowButtonMotionFcn',@mousemovepl3);
            
        else % if inside the plane: movement parallel to plane normal
            % it becomes a (2d) problem in the viewing plane. hitpoint shall be moved along the normal vector
            hitpr = rot*hitpoint'; 
            vecpr = rot*(hitpoint'+plotstruct.planedef{plane}.normal); 
            set(gcf,'WindowButtonMotionFcn',@mousemovepl1);
        end
        set(gcf,'pointer','hand');
        set(gcf,'WindowButtonUpFcn',@mouseuppl);
     end
     
        function moveok = moveplane(planedef) % helper function to do draw plane at new position
            % calculate new intersection with this plane with the axes box
           [planevectors, origin] = getplaneparameter(planedef.normal(:)', planedef.c);
           % intersection (like in doplot.m)
           [cutvertices, cutorder] = createmesh (corners, 1:8, planevectors', origin);
           if isempty(cutvertices), moveok=false; return; end % cannot move outside
           % transform projected ccordinates back to full system
           cutvertices = cutvertices * planevectors' + repmat(origin(:)',size(cutvertices,1),1);
           set(thisplaneh,'Vertices',cutvertices, 'Faces',cutorder);
           moveok = true;
        end
            
        function mousemovepl1(src,evnt) % To shift the plane along its normal
           showpos;
           cp = get(plotstruct.axeshandle,'CurrentPoint');
           cp2d = rot*cp(1,:)'; cp2d = cp2d(1:2)';
           newpr = hitpr + (cp2d-hitpr(1:2)')*(vecpr(1:2)-hitpr(1:2))/norm(vecpr(1:2)-hitpr(1:2))^2 * (vecpr-hitpr); % for the lambda, work only in 2d
           newp3d = rot \ newpr;
           newplanedef = plotstruct.planedef{plane};
           newplanedef.c = plotstruct.planedef{plane}.normal' * newp3d;
           if moveplane(newplanedef), plotstruct.planedef{plane} = newplanedef; end % move, and update if ok
        end
        
        function mousemovepl2(src,evnt) % To shift point along line
           showpos;
           cp = get(plotstruct.axeshandle,'CurrentPoint');
           cpr = rot*cp(1,:)';  mvdir=[0,0,0]; mvdir(mvdim)=1; mvdirr=rot*mvdir';
           mvpoint = [get(pointh(1),'xdata'), get(pointh(1),'ydata'), get(pointh(1),'zdata')];
           mvpointr = rot*mvpoint';
           newpr = mvpointr + ((cpr(1:2)-mvpointr(1:2))'*mvdirr(1:2))/norm(mvdirr(1:2)) * mvdirr;
           newp3d = rot \ newpr;
           mvpoint(mvdim) = newp3d(mvdim);
           newplanedef = plotstruct.planedef{plane};           
           [newplanedef.normal,newplanedef.c] = fithyperplane ([mvpoint; vertp(ncneighb,:)]);
           if ~moveplane(newplanedef), return;  end % move, and update if ok
           plotstruct.planedef{plane} = newplanedef;
           set(pointh(1),'xdata',mvpoint(1),'ydata',mvpoint(2),'zdata',mvpoint(3));
        end
        
         function mousemovepl3(src,evnt) % To shift edge
           showpos;
           cp = get(plotstruct.axeshandle,'CurrentPoint');           
           cpp = cp(1,:) + (cp(2,:)-cp(1,:))* (vertp(endp(1),cdim)-cp(1,cdim))/(cp(2,cdim)-cp(1,cdim)); % mousepoint on movement plane
           % cpp+ lambda*dvec is new line, intersect this with the box's face
           vertices = createmesh (axface2d, 1:4, dvec(setdiff(1:3,cdim)), cpp(setdiff(1:3,cdim))); % cut in 2d. * Maybe better use createintersection
           if size(vertices,1)~=2, return; end
           vertices = vertices*dvec + repmat(cpp,2,1); vertices(:,cdim) = vertp(endp(1),cdim);  % transform proj. back and set 3rd dim
           newplanedef = plotstruct.planedef{plane};           
           [newplanedef.normal,newplanedef.c] = fithyperplane ([vertices; vertp(nc,:)]);
           if ~moveplane(newplanedef), return;  end % move, and update if ok
           plotstruct.planedef{plane} = newplanedef;
           set(lineh,'xdata',vertices(:,1),'ydata',vertices(:,2),'zdata',vertices(:,3));
        end
   
        function mouseuppl(src,evnt)
           % Reset behaviour of figure
           set(src,'Pointer','arrow')
           set(src,'WindowButtonMotionFcn',@showpos)
           set(src,'WindowButtonUpFcn','')
           set(thisplaneh,'facecolor',oldcolor);
           delete(pointh); delete(lineh);
           set(plotstruct.axeshandle,'xlimmode',oldlimmode,'ylimmode',oldlimmode,'zlimmode',oldlimmode); axis(xyzlim);
           guidata(gcf, plotstruct);
           % ** do a new extracted plot
        end
              
   end %mouseclickonplane


%% Subfunction:  Prepare a scan

    function preparescan(src,evnt)
        plotstruct = guidata(gcf);
        plotstruct.scandef.type = 'interpolation';
        plotstruct.scandef.np = 20;  % Standard value for number of points, changed later
        plotstruct.scandef.xaxiscoord = 'AUTO';
        guidata(gcf, plotstruct);
        if threedim, createscanset; return; end % (in 3d mode, continue elsewhere)
        % Creates a second set of axes in the figure
        if isempty(plotstruct.scanaxeshandle)
            createscanaxes;
        end     
        % Set callback routine for mouse button 
        set(plotstruct.axeshandle,'ButtonDownFcn',@mouseclick1);
        set(get(plotstruct.axeshandle,'Children'),'ButtonDownFcn',@mouseclick1);
    end %preparescan

    function prepareintegration(src,evnt)
        plotstruct = guidata(gcf);
        plotstruct.scandef.type = 'integration';
        plotstruct.scandef.np = 20;  % Standard value for number of points, changed later
        plotstruct.scandef.xaxiscoord = 'AUTO';
        guidata(gcf, plotstruct);
        if threedim, createscanset; return; end % (in 3d mode, continue elsewhere)
        % Creates a second set of axes in the figure
        if isempty(plotstruct.scanaxeshandle)
            createscanaxes;
        end    
        % Set callback routine for mouse button 
        set(plotstruct.axeshandle,'ButtonDownFcn',@mouseclick2);
        set(get(plotstruct.axeshandle,'Children'),'ButtonDownFcn',@mouseclick2);
    end %prepareintegration

    function prepareprojection(src,evnt)
        % Creates a second set of axes in the figure
        plotstruct = guidata(gcf);
        plotstruct.scandef.type = 'projection';
        plotstruct.scandef.xaxiscoord = 'AUTO';
        guidata(gcf, plotstruct);
        if isempty(plotstruct.scanaxeshandle)
            createscanaxes;
        end    
        % Set callback routine for mouse button 
        set(plotstruct.axeshandle,'ButtonDownFcn',@mouseclick2);
        set(get(plotstruct.axeshandle,'Children'),'ButtonDownFcn',@mouseclick2);
    end %prepareprojection

    function createscanaxes
        % Resize the figure window conveniently...
            plotstruct = guidata(gcf);
            scrsz = get(0,'ScreenSize');    % Screen size
            figpos = get(gcf, 'Position');  % Figure window size
            figpos(3) = min(scrsz(3), 2*figpos(3));
            figpos(1) = min(figpos(1), scrsz(3)-figpos(3));
            figpos(4) = figpos(3) / 2;
            set(gcf, 'Position', figpos);
            set(plotstruct.axeshandle,'Outerposition',[0,0,.5,.95]');
            set(findobj('Tag','Titlebar2'),'Position',[0,.93,.5,.07]);
            set(plotstruct.poshandle, 'Position', [.3, 0, .195, .05]);
            drawnow;
            % Create new axes for scan plot
            plotstruct.scanaxeshandle = axes('OuterPosition',[.5,0,.5,.95]);
            menuentry = uimenu('Label','Scan');
                         uimenu(menuentry,    'Label', 'Scan definition',        'Callback', @inputscandef);
                         uimenu(menuentry,    'Label', 'Set number of points',   'Callback', callback_np);
            scanxaxis  = uimenu(menuentry,    'Label', 'Choose coordinate for x-axis');
                    if strcmpi(plotstruct.type, 'qeplane')
                            uimenu(scanxaxis, 'Label', 'Automatic',              'Callback', callback_xax_auto);
                            uimenu(scanxaxis, 'Label', '|Q|',                    'Callback', callback_xax_qm);
                            uimenu(scanxaxis, 'Label', 'En',                     'Callback', callback_xax_en);
                    end        
                    if strcmpi(plotstruct.type, 'qxy')        
                            uimenu(scanxaxis, 'Label', 'Automatic (Qx/Qy)',      'Callback', callback_xax_auto);
                            uimenu(scanxaxis, 'Label', '|Q|',                    'Callback', callback_xax_qm);
                            uimenu(scanxaxis, 'Label', 'Qx',                     'Callback', callback_xax_qx);
                            uimenu(scanxaxis, 'Label', 'Qy',                     'Callback', callback_xax_qy);
                            uimenu(scanxaxis, 'Label', 'H',                      'Callback', callback_xax_qh);
                            uimenu(scanxaxis, 'Label', 'K',                      'Callback', callback_xax_qk);
                            uimenu(scanxaxis, 'Label', 'L',                      'Callback', callback_xax_ql);
                    end
                         uimenu(menuentry,    'Label', 'Fit Gauss (nfit)',       'Callback', @gaussfit, 'Separator', 'on');
                         uimenu(menuentry,    'Label', 'Fit Linear (nfit)',      'Callback', @linearfit);
                         uimenu(menuentry,    'Label', 'Open nfit window',       'Callback', @generalnfit);
            scsavemenu = uimenu(menuentry,    'Label', 'Save scan data to file', 'Enable', 'off', 'Separator', 'on');  % Here, retain the handle to later enable it
                            uimenu(scsavemenu,'Label', 'Full format: (qx,qy,H,K,L,y,dy) + header', 'Callback', @savescan);
                            uimenu(scsavemenu,'Label', 'Simple format: as displayed (x, y, dy)  ', 'Callback', callback_scansave);
            scexportmen =uimenu(menuentry,    'Label', 'Export data to mfit window', 'Enable', 'off', 'Callback', @exportmfit);   
                         uimenu(menuentry,    'Label', 'Replot',                 'Callback', callcack_replotscan,'Separator', 'on' );
                         uimenu(menuentry,    'Label', 'Clear and reset window', 'Callback', @clearscan, 'Separator', 'on' );
            
%             if strcmpi(plotstruct.type,'qeplane')
%                 set(findobj('Label','Scan definition', 'Parent', menuentry), 'Enable', 'off');
%             end
                         
            % Save plotstruct
            guidata(gcf,plotstruct);
    end

%% Subfunction: Create a set of scans from multiple surfaces in 3D

    function createscanset
        % in 3d, the definition of a scan by mouse is not possible
        % Go directly to the numerical input.
        inputscandef([],[]); % (scan def. stored in plotstruct)
        
        plotstruct = guidata(gcf);
        
        if ~isfield(plotstruct.scandef,'xdat') || isempty(plotstruct.scandef.xdat) || ~isfield(plotstruct.scandef,'ydat') || isempty(plotstruct.scandef.ydat), return; end
        
        [planedef.normal, planedef.c] = fithyperplane([ plotstruct.scandef.xdat(1), plotstruct.scandef.ydat(1), 0; ...
                                                        plotstruct.scandef.xdat(2), plotstruct.scandef.ydat(2), 0; ...
                                                        plotstruct.scandef.xdat(1), plotstruct.scandef.ydat(1), 1] );        
%         % The separate sclices have to be identified
%         if ~isfield(plotstruct,'sectionlist') || isempty(plotstruct.sectionlist)
%             [~,section] = finddistinctsections(plotstruct.datalist.faces);
%             plotstruct.sectionlist = section;
%         end
%         section = plotstruct.sectionlist;
%         nsecs = max(section); 
        
        % Initialize new data structure for result
        newdata = plotstruct.datalist{1}; %.. on the basis of an existing, and delete contents
        newdata.coordlist = [];
        newdata.valuelist = [];
        newdata.dataname = ['Separate ', plotstruct.scandef.type, ' of multiple datasets'];
        newdata.raw = false;
        if isfield(newdata,'delaunaytri'), newdata = rmfield(newdata,'delaunaytri'); end
        newdata = rmfield(newdata,{'vertexlist','faces','monitorlist'});
        if strcmpi(plotstruct.type,'qxyz') 
            newdata.coordtype = 'qxy';
            % get vectors defining the new plane (cutplane)
            newdata.constants = [newdata.constants, 'QVERT'];
            newdata.QVERT = planedef.c;
            avec = getplaneparameter([planedef.normal'; 0,0,1],[planedef.c;0]); % horizontal vector in cutplane
            % ** try to make it point in right direction... ( to be generalized)
            if (abs(avec(1))>abs(avec(2)) && avec(1)<0) || (abs(avec(2))>abs(avec(1)) && avec(2)<0), avec=-avec; end
            bvec = cross(avec,planedef.normal);
            UB = UBmatrix(plotstruct.sampleinfo.lattice, plotstruct.sampleinfo.ax, plotstruct.sampleinfo.bx );
            newdata.sampleinfo.ax = makeinteger((UB\avec)',1e-4);
            newdata.sampleinfo.bx = makeinteger((UB\bvec)',1e-4);            
        elseif strcmpi(plotstruct.type,'qxqyen'), newdata.coordtype = 'qeplane'; 
        end
        
        [stdgrid] = getoption('stdgrid'); stdgrid = stdgrid.(upper(plotstruct.type));
        [ndata, ~, ~, pointinds, vertinds] = slicecounting(plotstruct.datalist);
        
        % Principle (simplified): use first two coordinates, ignore third. 
        % ** Have to generalize this!
        % ** Question: what to use as the x-axis?? 
        % ** (at least for 3d-Q cuts this is not always ok!)
        
        % For each of these slices
        tempplotstruct = plotstruct;
        tempplotstruct = rmfield(tempplotstruct,'datalist'); tempplotstruct.datalist{1} = plotstruct.datalist{1};
        tempplotstruct.type = 'QXY';
%         tempplotstruct.vertexlist = tempplotstruct.vertexlist(:,1:2);
        if isfield(tempplotstruct.datalist{1},'delaunaytri'), tempplotstruct.datalist{1} = rmfield(tempplotstruct.datalist{1},'delaunaytri'); end
        tempplotstruct.scandef.np = round(tempplotstruct.scandef.np / ndata);  % ?
        
        nonconstcount = 0;
        for nsl=1:ndata
            
            yval = plotstruct.coordlist(pointinds{nsl},3);
            if max(yval)-min(yval) > stdgrid(3), nonconstcount = nonconstcount+1; end
            yval = mean(yval);
            
            tempplotstruct.datalist{1}.faces = plotstruct.datalist{nsl}.faces;
            tempplotstruct.datalist{1}.coordlist = plotstruct.datalist{nsl}.coordlist;
            tempplotstruct.datalist{1}.valuelist = plotstruct.datalist{nsl}.valuelist;
            tempplotstruct.datalist{1}.monitorlist = plotstruct.datalist{nsl}.monitorlist;
            tempplotstruct.coordlist = plotstruct.coordlist(pointinds{nsl},1:2);
            tempplotstruct.datalist{1}.vertexlist = plotstruct.datalist{nsl}.vertexlist(:,1:2);
            tempplotstruct.vertexlist = plotstruct.vertexlist(vertinds{nsl},1:2);
                        
            [scandata, xlabeltext] = makescandata(tempplotstruct);
            
            validdata = isfinite(scandata.y);
            newdata.coordlist = [newdata.coordlist; scandata.x(validdata), yval*ones(sum(validdata),1)];
            newdata.valuelist = [newdata.valuelist; scandata.y(validdata), scandata.dy(validdata)];
            
        end
        
        if nonconstcount>0, fprintf('Warning: The third coordinate is not always constant. This is neglected in the calculation!\n'); end
                      
        % Add a this plane to plotstruct

        planedef.properties = {'FaceColor','r','FaceAlpha',.35,'Tag','cutplane','ButtonDownFcn',@mouseclickonplane}; %** use options.m for this?
        if ~isfield(plotstruct,'planedef'), plotstruct.planedef = {}; end
        plotstruct.planedef = [plotstruct.planedef, planedef];
        
        guidata(gcf, plotstruct);
        doplot(plotstruct);
        
         % ** Parameters to create the color plot?
        
        fcplot({newdata},newdata.coordtype);
        xlabel(xlabeltext);
        
    end

%% Subfunction: Analyze input of Scan definition

    function qxy = translateinput(text)
        
        function errormessage
            errordlg(['Format: In each line, give either two values in reciprocal Angstrom (for example: "qxy  -2.5  3.2") or three HKL (ex.: "hkl 2 0 0").' ...
                      'For q-en plots, give two values for |q| (A-1) and En (meV), "qen 1.5 2.0".']); 
        end
        
        qxy = [];
        qvwarning = false;
        for i = 1:length(text)
            [st] = regexpi(text{i},'[A-z]*');    % Get text  
            if numel(st)~=1 || length(text{i})<(st+2), errormessage; return; end
            switch upper(text{i}(st:(st+2)))
                case 'QEN', coord='qen';
                case 'QXY', coord='qxy';
                case 'HKL', coord='hkl'; 
                case 'THI', coord='thick'; if i<3, errormessage; return; end
                otherwise, errormessage; return
            end
            [st,en] = regexpi(text{i},'-?\d*\.?\d*');  % Get numbers
            switch coord
                case {'qxy', 'qen'} 
                    if numel(st)<2, errormessage; return; end  % Need 2 numbers
                    qxy{i} = [str2double(text{i}(st(1):en(1))), str2double(text{i}(st(2):en(2)))]; %#ok<AGROW>
                case 'hkl'
                    if numel(st)<3, errormessage; return; end  % Need 3 numbers
                    hkl = [str2double(text{i}(st(1):en(1))), str2double(text{i}(st(2):en(2))), str2double(text{i}(st(3):en(3)))];
                    UB = UBmatrix( plotstruct.sampleinfo.lattice, plotstruct.sampleinfo.ax, plotstruct.sampleinfo.bx );
                    [xh, yh, zh] = calcQS(1,0,0, UB);
                    [xk, yk, zk] = calcQS(0,1,0, UB);
                    [xl, yl, zl] = calcQS(0,0,1, UB);
                    qxy{i} = [hkl(1)*xh + hkl(2)*xk + hkl(3)*xl, hkl(1)*yh + hkl(2)*yk + hkl(3)*yl]; %#ok<AGROW>
                    % ** Check if this is really ok **
                    % ** {1} !
                    maxdeviate = getoption('maxdeviate');
                    if isfield(plotstruct.datalist{1},'QVERT') && abs((hkl(1)*zh + hkl(2)*zk + hkl(3)*zl)-plotstruct.datalist{1}.QVERT) > maxdeviate.QVERT  && ~qvwarning
                        warndlg('Attention: HKL entered are not consistent with vertical momentum transfer. Continue, but please check.','','modal');
                        qvwarning = true;
                    end
                case 'thick'  
                    if numel(st)<1, errormessage; return; end  % Need 1 number
                    vec = qxy{2}-qxy{1};
                    qxy{i} = [qxy{1}(1), qxy{1}(2)] + [-vec(2), vec(1)] * str2double(text{i}(st(1):en(1)))/2/sqrt(sum(vec.^2)); %#ok<AGROW>
            end
        end        
    end %translateinput
        

%% Callback: Numerical Input of Scan Definition

    function inputscandef(src,evnt)
        plotstruct = guidata(gcf);     
        plotopt  = getoption('plotopt');
        if plotopt.preferHKL   % use HKL instead Angstroms
            UB = UBmatrix( plotstruct.sampleinfo.lattice, plotstruct.sampleinfo.ax, plotstruct.sampleinfo.bx );
            try
                [H1,K1,L1] = calcHKL(plotstruct.scandef.xdat(1), plotstruct.scandef.ydat(1), plotstruct.datalist.QVERT, UB);
                [H2,K2,L2] = calcHKL(plotstruct.scandef.xdat(2), plotstruct.scandef.ydat(2), plotstruct.datalist.QVERT, UB);
                H1=round(H1*1e4)/1e4; K1=round(K1*1e4)/1e4; L1=round(L1*1e4)/1e4; H2=round(H2*1e4)/1e4; K2=round(K2*1e4)/1e4; L2=round(L2*1e4)/1e4;
                if any(strcmpi(plotstruct.scandef.type,{'Projection','Integration'}))
                    [H3,K3,L3] = calcHKL(plotstruct.scandef.xdat(1)+plotstruct.scandef.framevec(1), plotstruct.scandef.ydat(1)+plotstruct.scandef.framevec(2), plotstruct.datalist.QVERT, UB);
                    H3=round(H3*1e4)/1e4; K3=round(K3*1e4)/1e4; L3=round(L3*1e4)/1e4;
                end
            catch
                plotopt.preferHKL = false;
            end
        end
        if strcmpi(plotstruct.type, 'qxy'), stdcoord = 'qxy '; 
        elseif strcmpi(plotstruct.type, 'qeplane'), stdcoord = 'qen ';
        elseif strcmpi(plotstruct.type,'qxqyen'), stdcoord = 'qxy';
        end
        if strcmpi(plotstruct.scandef.type,'Interpolation')
            try
                if plotopt.preferHKL && ~strcmpi(plotstruct.type, 'qeplane')
                    defaulttext = {['hkl ' num2str(H1) ' ' num2str(K1) ' ' num2str(L1)],  ['hkl ' num2str(H2) ' ' num2str(K2) ' ' num2str(L2)] };
                else
                    defaulttext = {[stdcoord num2str(plotstruct.scandef.xdat(1)) '  ' num2str(plotstruct.scandef.ydat(1))],  ...
                                   [stdcoord num2str(plotstruct.scandef.xdat(2)) '  ' num2str(plotstruct.scandef.ydat(2))] };
                end
            catch
                defaulttext = {'', ''}; 
            end
            answ = inputdlg({'Start point', 'End point'}, 'Give scan definition', 1, defaulttext);
        elseif any(strcmpi(plotstruct.scandef.type,{'Projection','Integration'}))
            try
                if plotopt.preferHKL && ~strcmpi(plotstruct.type, 'qeplane')
                    defaulttext = {['hkl '  num2str(H1) ' ' num2str(K1) ' ' num2str(L1)], ...
                                   ['hkl '  num2str(H2) ' ' num2str(K2) ' ' num2str(L2)], ...
                                   ['hkl '  num2str(H3) ' ' num2str(K3) ' ' num2str(L3)]};

                else
                    defaulttext = {[stdcoord  num2str(plotstruct.scandef.xdat(1)) '  ' num2str(plotstruct.scandef.ydat(1))], ...
                                   [stdcoord  num2str(plotstruct.scandef.xdat(2)) '  ' num2str(plotstruct.scandef.ydat(2))], ...
                                   [stdcoord  num2str(plotstruct.scandef.xdat(1)+plotstruct.scandef.framevec(1)) '  ' num2str(plotstruct.scandef.ydat(1)+plotstruct.scandef.framevec(2))],};
                end
            catch
                defaulttext = {'', '', ''}; 
            end
            answ = inputdlg({'Start point', 'End point', 'Range: coord. of corner or thickness in axes units (ex.: "thick 0.5")'}, 'Give scan definition', 1, defaulttext);
        end
        
        plotstruct.scandef.xdat = [];
        plotstruct.scandef.ydat = [];
        
        if isempty(answ), return; end
        qxy = translateinput(answ);
        if isempty(qxy), return; end
%         plotstruct.scandef.xdat = [qxy{1}(1), qxy{2}(1)];
%         plotstruct.scandef.ydat = [qxy{1}(2), qxy{2}(2)];
        if any(strcmpi(plotstruct.scandef.type,{'Projection','Integration'})), framevec = qxy{3}-qxy{1}; else framevec = [0,0]; end
        if isfield(plotstruct,'qvert')
            xydat = makeliststruct([qxy{1}; qxy{2}; framevec],[],'coordtype','qxy','sampleinfo',plotstruct.sampleinfo,'QVERT',plotstruct.qvert);
        else
            xydat = makeliststruct([qxy{1}; qxy{2}; framevec],[],'coordtype','qxy','sampleinfo',plotstruct.sampleinfo,'QVERT',0);
        end
        % Transform to figure coordinate system:  ** make this more general.
        % if hkl:
        if strcmpi(plotstruct.type,'hklvectors'), xydat = coordtransform(xydat, 'hklvectors'); end
        if isempty(xydat), return; end
        plotstruct.scandef.xdat = xydat.coordlist(1:2,1)';
        plotstruct.scandef.ydat = xydat.coordlist(1:2,2)';    
        if any(strcmpi(plotstruct.scandef.type,{'Projection','Integration'}))
            plotstruct.scandef.framevec = xydat.coordlist(3,:);
        end
        
        [ndata, ~, ~, ~, vertinds] = slicecounting(plotstruct.datalist);
        plotstruct.scandef.np = 0;
        for nsl=1:ndata
            plotstruct.scandef.np = plotstruct.scandef.np + estimateNP([plotstruct.scandef.xdat', plotstruct.scandef.ydat'], plotstruct.vertexlist(vertinds{nsl},1:2), plotstruct.datalist{nsl}.faces);
        end
        
        if threedim, guidata(gcf, plotstruct); return; end % in 3D, stop here. Otherwise, continue with setting the lines for visualisation etc
        
        if strcmpi(plotstruct.scandef.type,'Interpolation')
            if isempty(plotstruct.scanlinehandle)
                axes(plotstruct.axeshandle);
                plotstruct.scanlinehandle = line('XData', plotstruct.scandef.xdat, 'YData', plotstruct.scandef.ydat, 'Marker','p','color','r');
            else
                set(plotstruct.scanlinehandle, 'XData', plotstruct.scandef.xdat, 'YData', plotstruct.scandef.ydat); drawnow 
            end
        elseif any(strcmpi(plotstruct.scandef.type,{'Projection','Integration'}))
            try %#ok<ALIGN>
            delete(plotstruct.scanlinehandle);
            delete(plotstruct.projectlinehandle);   
            delete(plotstruct.pp1handle);          
            delete(plotstruct.pp2handle);          
            delete(plotstruct.projectframehandle); 
            catch end
            axes(plotstruct.axeshandle);
            plotstruct.scanlinehandle=[];
            plotstruct.projectlinehandle = line('XData', plotstruct.scandef.xdat, 'YData', plotstruct.scandef.ydat,'Marker','none','color','r','ButtonDownFcn', @mouseclick_projectline);
            plotstruct.projectframehandle = line('XData', [plotstruct.scandef.xdat(1)+framevec(1), plotstruct.scandef.xdat(2)+framevec(1), plotstruct.scandef.xdat(2)-framevec(1), ...
                                                           plotstruct.scandef.xdat(1)-framevec(1), plotstruct.scandef.xdat(1)+framevec(1)], ...
                                                 'YData', [plotstruct.scandef.ydat(1)+framevec(2), plotstruct.scandef.ydat(2)+framevec(2), plotstruct.scandef.ydat(2)-framevec(2), ...
                                                           plotstruct.scandef.ydat(1)-framevec(2), plotstruct.scandef.ydat(1)+framevec(2)],'Marker','none','color','r');
            plotstruct.pp1handle = line('XData',plotstruct.scandef.xdat(2),'YData',plotstruct.scandef.ydat(2),'Marker','o','MarkerFaceColor','r','color','r', 'ButtonDownFcn', @mouseclick_pp1);
            plotstruct.pp2handle = line('XData',plotstruct.scandef.xdat(1)+framevec(1),'YData',plotstruct.scandef.ydat(1)+framevec(2),'Marker','o','MarkerFaceColor','r','color','r', 'ButtonDownFcn', @mouseclick_pp2);
            
        end       
            
        guidata(gcf, plotstruct);
        doscanplot(plotstruct);
        set(scsavemenu, 'Enable', 'on');
        set(scexportmen, 'Enable', 'on');

    end %inputscandef




%% Callback: Delete scan axes and reset

    function clearscan(src,evnt)
        plotstruct = guidata(gcf);
        if isempty(plotstruct.scanaxeshandle), return; end
        delete(plotstruct.scanaxeshandle);
        plotstruct.scanaxeshandle = [];
        delete(findobj('Type','uimenu','Label','Scan'));
        try  %#ok<ALIGN> % Delete, if present
            delete(plotstruct.scanlinehandle);     
            delete(plotstruct.projectlinehandle);  
            delete(plotstruct.pp1handle);           
            delete(plotstruct.pp2handle);           
            delete(plotstruct.projectframehandle); 
        catch end
        plotstruct.scanlinehandle=[]; plotstruct.projectlinehandle=[]; plotstruct.pp1handle=[]; plotstruct.pp2handle=[]; plotstruct.projectframehandle=[];
        figpos = get(gcf, 'Position');  % Figure window size
        figpos(3) = 1/2 * figpos(3);
        set(gcf, 'Position', figpos);
        set(plotstruct.axeshandle,'Outerposition',[0,0,1,.95]');
        set(plotstruct.poshandle, 'Position', [.6, 0, .39, .05]);
        drawnow;
        
        % Remove callback routine for mouse button 
        set(plotstruct.axeshandle,'ButtonDownFcn',[]);
        set(get(plotstruct.axeshandle,'Children'),'ButtonDownFcn',[]);
        
        guidata(gcf, plotstruct);
    end

%% Callback: Delete a point by clicking

    function start_clickdelete(src,evnt)  % Called by menu, sets the mouse callback routine
        plotstruct = guidata(gcf); 
        set(gcf,'Pointer','crosshair')
        set(plotstruct.axeshandle,'ButtonDownFcn',@deleteclick); 
        set(get(plotstruct.axeshandle,'Children'),'ButtonDownFcn',@deleteclick) ;
        % Ensure existence of delaunay
        for nsl=1:length(plotstruct.datalist)
            if ~isfield(plotstruct.datalist{nsl},'delaunaytri')
                stdratio = getoption('stdratio'); stdratio = stdratio.(upper(plotstruct.datalist{nsl}.coordtype));
                % ** to be generalized for 3d!
                try
                plotstruct.datalist{nsl}.delaunaytri = delaunay(stdratio(1) * plotstruct.datalist{nsl}.coordlist(:,1), stdratio(2) * plotstruct.datalist{nsl}.coordlist(:,2));
                catch
                plotstruct.datalist{nsl}.delaunaytri = delaunay(stdratio(1) * plotstruct.datalist{nsl}.coordlist(:,1), stdratio(2) * plotstruct.datalist{nsl}.coordlist(:,2),{'Qt','Qbb','Qc','Qz'});
                end
            end
        end
        guidata(gcf,plotstruct);
    end %start_clickdelete


    function deleteclick(src,evnt)
        plotstruct = guidata(gcf);
        
        % If right mouse button pressed, quit delete mode
        if strcmpi(get(gcf,'SelectionType'),'alt')  
            set(plotstruct.axeshandle,'ButtonDownFcn',[]); 
            set(get(plotstruct.axeshandle,'Children'),'ButtonDownFcn',[]) ;
            set(gcf,'Pointer','arrow');
            % Eventually set callback routines for movements
            set(plotstruct.projectlinehandle, 'ButtonDownFcn', @mouseclick_projectline);
            set(plotstruct.pp1handle, 'ButtonDownFcn', @mouseclick_pp1);
            set(plotstruct.pp2handle, 'ButtonDownFcn', @mouseclick_pp2);
            set(plotstruct.scanlinehandle, 'ButtonDownFcn', @mouseclick_line);  
            set(findobj(get(plotstruct.axeshandle,'children'),'Tag','cutplane'), 'ButtonDownFcn', @mouseclickonplane);
            return
        end
        
        % Otherwise, delete the clicked face
%         cp = get(plotstruct.axeshandle,'CurrentPoint');
%         nface = findface(cp(1,1), cp(1,2), plotstruct.coordlist, plotstruct.datalist.faces, plotstruct.vertexlist, plotstruct.datalist.delaunaytri);
        [~,nslice,nfaceonslice] = findface3d(plotstruct.axeshandle, plotstruct.coordlist, plotstruct.vertexlist, plotstruct.datalist); % use new findface3d instead findface
        
        plotstruct = deletepoint(plotstruct,nfaceonslice,nslice);
        
        guidata(gcf, plotstruct);
        doplot(plotstruct, 'keepview');
       
    end %deleteclick

%% Callback: Select a region to delete
    function start_deleteregion(src,evnt)
        plotstruct = guidata(gcf);
        set(gcf,'Pointer','crosshair')
        set(plotstruct.axeshandle,'ButtonDownFcn',@selectregion); 
        set(get(plotstruct.axeshandle,'Children'),'ButtonDownFcn',@selectregion) ;

        
        function selectregion(src,evnt)
            plotstruct = guidata(gcf);
            
            % If right mouse button pressed, quit delete mode
            if strcmpi(get(gcf,'SelectionType'),'alt')  
                set(plotstruct.axeshandle,'ButtonDownFcn',[]); 
                set(get(plotstruct.axeshandle,'Children'),'ButtonDownFcn',[]) ;
                set(gcf,'Pointer','arrow');
                % Eventually set callback routines for movements
                set(plotstruct.projectlinehandle, 'ButtonDownFcn', @mouseclick_projectline);
                set(plotstruct.pp1handle, 'ButtonDownFcn', @mouseclick_pp1);
                set(plotstruct.pp2handle, 'ButtonDownFcn', @mouseclick_pp2);
                set(plotstruct.scanlinehandle, 'ButtonDownFcn', @mouseclick_line);
                set(findobj(get(plotstruct.axeshandle,'children'),'Tag','cutplane'), 'ButtonDownFcn', @mouseclickonplane);   
                return
            end
            
            % Otherwise, go on
            % identify which slice clicked and work only on this one in the follwing
            [~,nslice] = findface3d(plotstruct.axeshandle, plotstruct.coordlist, plotstruct.vertexlist, plotstruct.datalist);
            if ~isfinite(nslice)
                if length(plotstruct.datalist)==1, nslice=1; else return; end
            end
            
            if ~isfield(plotstruct.datalist{nslice},'delaunaytri')
                % ** generalize this for 3d!
                stdratio = getoption('stdratio'); stdratio = stdratio.(upper(plotstruct.datalist{nslice}.coordtype));
                try
                plotstruct.datalist.delaunaytri{nslice} = delaunay(stdratio(1) * plotstruct.datalist{nslice}.coordlist(:,1), stdratio(2) * plotstruct.datalist{nslice}.coordlist(:,2));
                catch
                    fprintf('Abort - delaunay not yet implemented here.\n'); % **
                end
            end
            % Fraw rectangle
            % first corner of rectangle
            cp1 = get(plotstruct.axeshandle,'CurrentPoint');
            % draw and wait for button release
            rbbox;   % (Matlab command to drag rectangle)
            cp2 = get(plotstruct.axeshandle,'CurrentPoint');    % button up detected
            
            rot = rotationtoviewingplane(plotstruct.axeshandle);  % Get the rotation matrix to transform coords to viewing system
            cp1 = rot * cp1'; 
            cp2 = rot * cp2';
            corners = [cp1(1:2); cp2(1:2)];
            % extract polygon definition
            xp = [min(corners(:,1)), max(corners(:,1)), max(corners(:,1)), min(corners(:,1)), min(corners(:,1))];
            yp = [max(corners(:,2)), max(corners(:,2)), min(corners(:,2)), min(corners(:,2)), max(corners(:,2))];
            
            % find points inside rectangle
            [~,~,~,pointinds] = slicecounting(plotstruct.datalist);
            coordpoints = plotstruct.coordlist(pointinds{nslice},:); if size(coordpoints,2)==2, coordpoints(:,3)=0; end
            coordpoints = rot * coordpoints';
            incoord = inpolygon(coordpoints(1,:),coordpoints(2,:), xp, yp);
            
            % delete them
            plotstruct = deletepoint(plotstruct, find(incoord(:)), nslice);
            
            guidata(gcf, plotstruct);
            doplot(plotstruct, 'keepview');

        end
        
    end % select_deleteregion

%%  Callback: Change coordinate axes
    function changeaxes(src,evnt) 
        clearscan
        plotstruct = guidata(gcf);
 
        if strfind(upper(get(src,'Label')),'ANGLES'), changeto = 'ANGLES';
            % ** {1} !
        elseif strfind(upper(get(src,'Label')),'QX'), changeto = 'QXY'; if any(strcmpi(plotstruct.datalist{1}.coordtype,{'AnglesEnergy','QxQyEn'})), changeto = 'QXQYEN'; end
        elseif strfind(upper(get(src,'Label')),'LATTICE'), changeto = 'hklvectors';
        end
            
        plotstruct.type = changeto;     % ** {1} !
        [ctr, opt] = coordtransform(plotstruct.datalist{1},changeto);
        if isempty(ctr), fprintf('No plot - conversion into desired coordinate axes not successful.\n'); plotstruct = guidata(gcf);  return; end  
        plotstruct.vertexlist = ctr.vertexlist;
        plotstruct.coordlist = ctr.coordlist;
        if isfield(plotstruct,'linedef'), plotstruct.linedef = {};  end % Delete evtl. lines (** could also be transformed...)
        plotstruct.axesnames = readinput('axesnames',opt);
        if strcmpi(plotstruct.type,'hklvectors'), plotstruct.basisvectors = readinput('basisvectors',opt); elseif isfield(plotstruct,'basisvectors'), plotstruct=rmfield(plotstruct,'basisvectors'); end 
        guidata(gcf,plotstruct);
        doplot(plotstruct);
    end


%% Subfunction: Save scan to file

    function savescan(src,evnt)
        plotstruct = guidata(gcf);
        if strcmpi(plotstruct.scandef.type,'projection') || ~any(strcmpi(plotstruct.type,{'qxy','angles','hklvectors'}))
            msgbox('Long format only possible for interpolation and integration of constant energy data. Please use short format.','Error','error','modal'); return; 
        end
        
        [hbar, mass_n, meVJ, normalizeto, normval] = getoption('hbar', 'mass_n', 'meVJ', 'normalizeto', 'normval');
        % get filename from dialog box
        [file,path] = uiputfile('*.*','Save Flatcone scan data (Long format)'); 
        if file==0, return; end
         % get qx,qy and calculate HKL values
        startpoint = [plotstruct.scandef.xdat(1), plotstruct.scandef.ydat(1)];
        endpoint   = [plotstruct.scandef.xdat(2), plotstruct.scandef.ydat(2)];
        np         =  plotstruct.scandef.np;
        if isfield(plotstruct.scandef,'nslice'), sslice = plotstruct.scandef.nslice; else sslice = 1; end
        scanpath = repmat(startpoint,np,1) + (0:(np-1))'/np  * (endpoint-startpoint);
        UB = UBmatrix( plotstruct.sampleinfo.lattice, plotstruct.sampleinfo.ax, plotstruct.sampleinfo.bx);
        [H, K, L] = calcHKL( scanpath(:,1), scanpath(:,2), plotstruct.datalist{sslice}.QVERT, UB );      
        
        sfile = fopen([path,file], 'w');
        
        if sfile==-1, msgbox('Error on opening the file. Output could not be written.','Error','error','modal'); return; end
        
        % Write header
        fprintf(sfile, '%s %s\n', '# Scan from Flatcone data ', plotstruct.scandef.type);
        fprintf(sfile, '%s\n',['# ' plotstruct.datalist{sslice}.expname ', Data: ' plotstruct.datalist{sslice}.dataname ]);
        fprintf(sfile, '%s\n',['# Data normalized to ' normalizeto ' = ' num2str(normval)]); 
        fprintf(sfile, '%s\n',['# Vertical momentum transfer ' num2str(plotstruct.datalist{sslice}.QVERT) ' A^-1']);
        fprintf(sfile, '%s\n',['# Energy transfer ' num2str((plotstruct.datalist{sslice}.KI^2-plotstruct.datalist{sslice}.KF^2)*1E20*hbar^2/2/mass_n*meVJ) ' meV']);
        if any(strcmpi(plotstruct.datalist{sslice}.constants, 'TEMP')),  fprintf(sfile, '%s\n',['# Temperature ' num2str(plotstruct.datalist{sslice}.TEMP,4) ' K']);
        end
        fprintf(sfile, '%s ',['# Scan data obtained by ' plotstruct.scandef.type ' of original data.']);
        if any(strcmpi(plotstruct.scandef.type,{'integration','projection'}))
            fprintf(sfile,'%s ',['Width of ' plotstruct.scandef.type ' region in axes units ']);
            if strcmpi(plotstruct.type,'qxy'), fprintf(sfile,'(A-1)'); 
            elseif strcmpi(plotstruct.type,'angles'), fprintf(sfile,'(degrees)'); 
            elseif strcmpi(plotstruct.type,'hklvectors'), fprintf(sfile,'(vectors in rec. space)'); 
            else fprintf(sfile,'(mixed)'); 
            end
            fv = plotstruct.scandef.framevec; s=[plotstruct.scandef.xdat(2)-plotstruct.scandef.xdat(1), plotstruct.scandef.ydat(2)-plotstruct.scandef.ydat(1)];
            fprintf(sfile,' %6.3f', 2 * sqrt(fv*fv') * sqrt(1 - (fv*s')^2/(fv*fv')/(s*s')) );
            if abs(fv*s')>1e-5, fprintf(sfile,' (non-orthogonal frame)'); end
        end
            

        fprintf(sfile, '\n# \n%s \n', '# Qx(A-1) Qy(A-1)  H      K      L       signal      error');
        % Now write the data
        for i=1:np
            fprintf(sfile, ' %6.3f  %6.3f  %6.3f %6.3f %6.3f %10.2f %10.2f \n', scanpath(i,1), scanpath(i,2), H(i), K(i), L(i), plotstruct.scandata.y(i), plotstruct.scandata.dy(i));
        end
        
        fclose(sfile);
        
        guidata(gcf, plotstruct);
    end

%% Subfunction: Delete a point from plotstruct

    function plotstruct = deletepoint(plotstruct,nface,nslice)
        % helper function called by callback routines to delete one or
        % several points from the plotstruct
        
        if ~any(~isnan(nface)) || isempty(nface) ,  return ; end
        if nargin<3 || isempty(nslice), nslice=1; end;
        
%         % Determine the faces that share a border (i.e. two! vertices) with this one
%         vert = plotstruct.datalist.faces(nface, isfinite(plotstruct.datalist.faces(nface, :)));  % Vertices belonging to nface
%         eqv = (plotstruct.datalist.faces == vert(1));
%         for i=2:numel(vert)
%             eqv = eqv | (plotstruct.datalist.faces == vert(i)); 
%         end
%         neighbors = find (sum(eqv,2) >= 2);
%         neighbors = setdiff(neighbors, nface);

        [stdcell,stdratio,stdbindist] = getoption('stdcell','stdratio','stdbindist');
        stdratio = stdratio.(upper(plotstruct.datalist{nslice}.coordtype));
        maxsize = stdcell.(upper(plotstruct.datalist{nslice}.coordtype));
        stdbindist = stdbindist.(upper(plotstruct.datalist{nslice}.coordtype));
        % Check if can be treated in 2D mode (constant dims if d>2...) ** Generalize this ?
        vardims = (max(plotstruct.datalist{nslice}.coordlist,[],1)- min(plotstruct.datalist{nslice}.coordlist,[],1)) >= stdbindist(:)';
        constdim =  find(~vardims,1,'first');
        if sum(vardims)~=2, fprintf('Deleting not yet implemented for this type of 3d slice\n'); return; end
        
        [~, oldnpoints, oldnverts] = slicecounting(plotstruct.datalist);

        neighbors = findneighbors(nface, plotstruct.datalist{nslice}.coordlist, plotstruct.datalist{nslice}.delaunaytri);
        
        
       
        
        if numel(neighbors)>0  % if no neighbors, no other cells need to be recomputed
        
%             % Now determine also second neighbors (neighbors of neighbors)
%             ne_vert = plotstruct.datalist.faces(neighbors, :);
%             ne_vert = unique(ne_vert(isfinite(ne_vert)))';
%             eqv = (plotstruct.datalist.faces == ne_vert(1));
%             for i=2:numel(ne_vert)
%                 eqv = eqv | (plotstruct.datalist.faces == ne_vert(i)); 
%             end
%             ne_neighbors = find (sum(eqv,2) >= 2);
            
            ne_neighbors = findneighbors(neighbors, plotstruct.datalist{nslice}.coordlist, plotstruct.datalist{nslice}.delaunaytri);
            
            ne_neighbors = setdiff(ne_neighbors, [neighbors; nface]); 
            neighbors    = setdiff(neighbors, nface);

            % Now, to do the local reconstruction of the Voronoi diagram,
            % recompute it for neighbors and ne_neighbors. The cells for
            % neighbors will be correct, the ones for ne_neighbors not - they
            % serve to bound correctly the 'neighbors' in all directions.
            % Important: do the calculation in the coordinate space of
            % datalist!! (not in that of plot axes!!) Transform afterwards.
            
            [repairvertexlist,repairfaces,repairdelaunay] = ... 
                makevoronoi( plotstruct.datalist{nslice}.coordlist([neighbors;ne_neighbors],vardims), maxsize(vardims), stdratio(vardims), 'calctri');
            % evtl. add 3rd dimension
            if any(~vardims), 
                repairvertexlist = [repairvertexlist(:,1:constdim-1), ones(size(repairvertexlist,1),1)*mean(plotstruct.datalist{nslice}.coordlist(:,constdim)), repairvertexlist(:,constdim:end)]; end

            % Keep only neighboring faces, and only triangles between neighboring points
            repairfaces = repairfaces(1:numel(neighbors),:);
            repairdelaunay = repairdelaunay(all(ismember(repairdelaunay, 1:numel(neighbors)), 2), :);
                        
%             % Identify which are the new vertices appearing by the Voronoi reconstrucion.
%             % To be used below to repair Delaunay triangulation.
%             oldfaces = plotstruct.datalist.faces(neighbors,:); 
%             oldvertnum = unique(oldfaces(:)); oldvertnum = oldvertnum(isfinite(oldvertnum));
%             oldvertices = plotstruct.datalist.vertexlist(oldvertnum,:);
%             newvertnum = unique(repairfaces(:)); newvertnum = newvertnum(isfinite(newvertnum));
%             newvertices = repairvertexlist(newvertnum,:);
%             [diffvertices,diffvertind] = setdiff(round(newvertices*1E6), round(oldvertices*1E6), 'rows');  % ** Limit precision
%             diffvertices = newvertices(diffvertind,:);
%             % diffvertices now contains the newly appeared nodes.
            
            
            %Increase indices of vertices, before repairvertexlist will be appended to vertexlist
            repairfaces = repairfaces + size(plotstruct.datalist{nslice}.vertexlist,1); 
            
            % Replace changed faces (neighbors)
            repairfaces(:, (end+1):size(plotstruct.datalist{nslice}.faces,2) ) = NaN; % fill NaNs to overwrite existing entries
            plotstruct.datalist{nslice}.faces(neighbors, 1:size(repairfaces,2)) = repairfaces;
            plotstruct.datalist{nslice}.faces(plotstruct.datalist{nslice}.faces==0) = NaN;
            % Append new vertices
            plotstruct.datalist{nslice}.vertexlist = [plotstruct.datalist{nslice}.vertexlist; repairvertexlist]; 
            % Optimize patch (necessary, as in case of non-unique vertices neighbors may not be recognized)
            [plotstruct.datalist{nslice}.vertexlist, plotstruct.datalist{nslice}.faces] = optimizepatch(plotstruct.datalist{nslice}.vertexlist, plotstruct.datalist{nslice}.faces);
            
            % Transform this list
            ctr = coordtransform(plotstruct.datalist{nslice},plotstruct.type);
            % Insert at the right place in plotstruct.vertexlist
            plotstruct.vertexlist = [plotstruct.vertexlist(1:sum(oldnverts(1:nslice-1)),:); ctr.vertexlist; plotstruct.vertexlist(sum(oldnverts(1:nslice))+1:end,:)];
            
            
%              plotstruct.vertexlist = plotstruct.vertexlist(vassign,:);
        
        
                
            % Concept to repair Delaunay triangulation:
            % Delete all triangles that contain a deleted point.
            % Then use the information from the Voronoi cells:
            % Every corner of a Voronoi cell corresponds to exactly one center
            % of a triangle-enclosing circle ("Umkreis") (may contain several triangles, though). Therefore, 
            % the NEW Voronoi cells (has been done above), and add triangles for each cornerpoint.

            % Delete all triangles that contain a deleted point:
            plotstruct.datalist{nslice}.delaunaytri = plotstruct.datalist{nslice}.delaunaytri( ~any(ismember(plotstruct.datalist{nslice}.delaunaytri, nface),2), :);

            % Add new ones:
            neighbors = neighbors(:)';
            repairdelaunay = neighbors(repairdelaunay);   % convert delaunaytri to indices into real coordlist
            plotstruct.datalist{nslice}.delaunaytri = [plotstruct.datalist{nslice}.delaunaytri; repairdelaunay];
            
%             addtri = delaunayfromvoronoi( repairfaces, plotstruct.datalist.coordlist(neighbors,:) );
%             plotstruct.datalist.delaunaytri = [plotstruct.datalist.delaunaytri; addtri];
        
        end   
           
        
        
        
        
%         % ** ist das folgende wirklich richtig? Und kann man es effizienter machen?
%         % Oder einfach die Triangulation löschen?
%         nface = sort(nface,'descend');
%         for nf = 1:numel(nface)
%             % Repair Delaunay triangulation...!! (Local reconstruction!)
%             ind = sum(plotstruct.datalist.delaunaytri == nface(nf), 2) > 0; % lines (triangles) containig nface (logical index)
%             % Find all adjacent points
%             adj = plotstruct.datalist.delaunaytri(ind,:);
%             adj = unique( setdiff(adj(:)', nface(nf)) );
%             calcnewtri = false; 
%             
%             if numel(adj)> 2, try  %#ok<ALIGN>
%                 % Compute new triangulation for adjacent points
%                 newtri = delaunay(stdratio(1) * plotstruct.datalist.coordlist(adj,1), stdratio(2) * plotstruct.datalist.coordlist(adj,2));
%                 newtri = adj(newtri);  % (Attention, newtri was indices into limited list before)
%                 % Remove all triangles that contain the deleted point (keep
%                 % the lines that do not contain nface), and
%                 % Add new triangles
%                 plotstruct.datalist.delaunaytri = [plotstruct.datalist.delaunaytri(ind==0, :); newtri];
%                 
%             catch % In rare cases this does not work (insufficient number of adjacent points, bad configuration...)
%                   % In that case, recompute the whole delaunay triangulation
%                 plotstruct.datalist.delaunaytri = delaunay(stdratio(1) * plotstruct.datalist.coordlist([1:(nface(nf)-1),(nface(nf)+1):end],1), ...
%                                                            stdratio(2) * plotstruct.datalist.coordlist([1:(nface(nf)-1),(nface(nf)+1):end],2));
%                 calcnewtri = true;
%             end, end
%             
%             if ~calcnewtri, plotstruct.datalist.delaunaytri(plotstruct.datalist.delaunaytri > nface(nf)) = plotstruct.datalist.delaunaytri(plotstruct.datalist.delaunaytri > nface(nf)) - 1; end
%          
%         end
        
        lineind = setdiff(1:size(plotstruct.datalist{nslice}.coordlist,1), nface);
        
        linenum = nan(size(plotstruct.datalist{nslice}.coordlist,1),1); linenum(lineind) = 1:numel(lineind);
        plotstruct.datalist{nslice}.delaunaytri     = linenum(plotstruct.datalist{nslice}.delaunaytri);  % Renumber delaunaytri

        % Delete nface from everywhere...  // before, store the coordinates in deletelist 
        
        if isfield(plotstruct,'deletelist'), plotstruct.deletelist = [plotstruct.deletelist; plotstruct.datalist{nslice}.coordlist(nface,:)];
        else plotstruct.deletelist = plotstruct.datalist{nslice}.coordlist(nface,:); end
                
        plotstruct.datalist{nslice}.coordlist       = plotstruct.datalist{nslice}.coordlist(lineind,:);
        plotstruct.datalist{nslice}.valuelist       = plotstruct.datalist{nslice}.valuelist(lineind,:);
        if isfield(plotstruct.datalist{nslice},'monitorlist')
            plotstruct.datalist{nslice}.monitorlist = plotstruct.datalist{nslice}.monitorlist(lineind,:);
        end
        plotstruct.datalist{nslice}.faces           = plotstruct.datalist{nslice}.faces(lineind,:);
        
        plotstruct.coordlist                        = [plotstruct.coordlist(1:sum(oldnpoints(1:nslice-1)),:); ...
                                                       plotstruct.coordlist(sum(oldnpoints(1:nslice-1))+lineind,:); ...
                                                       plotstruct.coordlist(sum(oldnpoints(1:nslice))+1:end,:)];
        
        % If an interpolation had been explicitly calculated, it is now obsolete. Delete it.
        % ** Better delete only that for the selected slice
        if isfield(plotstruct,'interpolatedimage'), plotstruct = rmfield(plotstruct,'interpolatedimage'); end
        
        % Delete other things for that slice
        if isfield(plotstruct,'sliceedges') && length(plotstruct.sliceedges)>=nslice, plotstruct.sliceedges{nslice}=[]; end
        
        
    end %deletepoint

%% Callback: Set interpolation options
    function setinterpolationoptions(scr,evnt)
        plotstruct = guidata(gcf);
        % Construct default text from current values
        switch plotstruct.interpolationtype
            case 'pcolorsmooth', answ{1}='2';
            case 'pcolorfacet', answ{1}='1';
            case 'patchinterp', answ{1}='3';
        end
        answ{2}=plotstruct.interpolationalgorithm;
        answ{3}=num2str(plotstruct.interpolationgrid);
        answ{4}=num2str(plotstruct.interpolationlimit);
        answ{5}=plotstruct.interpolationsystem;
        % Call dialogbox
        answ = inputdlg({'Interpolation type: (1):Grid-pixelized, (2):Grid-smoothed, (3):Triangulation-smoothed', ...
                         'Interpolation algorithm: ''linear'',''cubic'',''natural'',''nearest''', ...
                         'Grid size (pixel per dim.) (attention, large values time consuming)', ...
                         'Relative distance limit for interpolation', ...
                         'Use coordinate system of original ''data'' or ''plot'':'}, ...
                        'Set options for interpolation', 1, answ);
        if isempty(answ), return; end
        readerror = false;
        switch strtrim(answ{1})
            case '1', plotstruct.interpolationtype = 'pcolorfacet';
            case '2', plotstruct.interpolationtype = 'pcolorsmooth';
            case '3', plotstruct.interpolationtype = 'patchinterp';
            otherwise readerror=true;
        end
        if ~any(strcmpi(strtrim(answ{2}),{'linear','cubic','natural','nearest'})), readerror = true;
        else plotstruct.interpolationalgorithm = lower(strtrim(answ{2})); end
        plotstruct.interpolationgrid = round(str2double(answ{3})); 
        if ~isfinite(plotstruct.interpolationgrid), readerror = true; end
        plotstruct.interpolationlimit = str2double(answ{4}); 
        if ~isfinite(plotstruct.interpolationlimit), readerror = true; end
        if ~any(strcmpi(strtrim(answ{5}),{'data','plot'})), readerror = true;
        else plotstruct.interpolationsystem = lower(strtrim(answ{5})); end
            
        if readerror, errordlg('Error: Input could not be evaluated.','Interpolation options'); return; end    
        
        if isfield(plotstruct,'interpolatedimage'), plotstruct = rmfield(plotstruct,'interpolatedimage'); end
        guidata(gcf, plotstruct);
        doplot(plotstruct);
    end



%% Callback: Choose Powder line to plot
    function show_powder(src,evnt)   % Choose and plot powder lines
        plotstruct = guidata(gcf);
        [dvals,proc] = powderlineselection;     % Open powder line selection window
        if numel(dvals)==1 && dvals==0, return; end   % Cancel was pressed, do no changes
        % Create field if not yet existing
        if ~isfield(plotstruct,'linedef'), plotstruct.linedef =  {}; end
        for i = length(plotstruct.linedef):-1:1
            if strcmpi(plotstruct.linedef{i}.label,'powder'), plotstruct.linedef = plotstruct.linedef([1:(i-1),(i+1):end]); end
        end % Clear powder lines
        % Set some plot style
        thislinedef.properties = {'edgecolor','k','linestyle','-'};
        thislinedef.label = 'powder';
        for nsl = 1:length(plotstruct.datalist) % loop over slices
            thislinedef.assigntoslice = nsl;
            for nl = 1:numel(dvals)
                for pr = find(proc)
                    % obtain intersection points of line with cell boundaries
                    dl = plotstruct.datalist{nsl}; % (use as "frame" structure for line's coords)
                    [dl.coordlist,thislinedef.connection] = calcpowderline(plotstruct.datalist{nsl}, dvals(nl), pr); 
                    if isempty(dl.coordlist), continue; end
                    % transform into plot coordinates
                    dl = coordtransform (dl, plotstruct.type); 
                    thislinedef.points = dl.coordlist;
                    % Add to list of lines
                    plotstruct.linedef = [plotstruct.linedef, thislinedef];
                end
            end
        end
        guidata(gcf, plotstruct);
        doplot(plotstruct);
    end

%% Callback: Draw intersections with HKL planes on surfaces
    function callback_HKLintersect(src,evnt)
        plotstruct = guidata(gcf);
        if ~isfield(plotstruct,'sampleinfo'); return; end %(only possible if hkl can be calculated)
        % get the Q-coordinates of the corners of the visible range
        xl = xlim(plotstruct.axeshandle); yl = ylim(plotstruct.axeshandle); zl = zlim(plotstruct.axeshandle);
        UB = UBmatrix( plotstruct.sampleinfo.lattice, plotstruct.sampleinfo.ax, plotstruct.sampleinfo.bx);
        switch upper(plotstruct.type)
            case {'QXY','QXQYEN'}
                extremepoints = [xl(1),yl(1); xl(2),yl(1); xl(1),yl(2); xl(2),yl(2) ];
                if numel(unique(round(100*mergelist(plotstruct.datalist,'QVERT'))/100))==1 %#ok<ALIGN> % check if all QVERT the same
                    extremepoints(:,3) = plotstruct.datalist{1}.QVERT;
                else return; end % not possible
                [H, K, L] = calcHKL( extremepoints(:,1), extremepoints(:,2), extremepoints(:,3), UB );
            case 'QXYZ'
                extremepoints = [xl(1),yl(1); xl(2),yl(1); xl(1),yl(2); xl(2),yl(2) ];
                extremepoints = [ extremepoints, [zl(1);zl(1);zl(1);zl(1)]; extremepoints, [zl(2);zl(2);zl(2);zl(2)] ];
                [H, K, L] = calcHKL( extremepoints(:,1), extremepoints(:,2), extremepoints(:,3), UB );
            case 'HKLVECTORS'
                H = plotstruct.basisvectors.origin(1) + [xl(1);xl(2);xl(1);xl(2)]*plotstruct.basisvectors.vector1(1) + [yl(1);yl(1);yl(2);yl(2)]*plotstruct.basisvectors.vector2(1);
                K = plotstruct.basisvectors.origin(2) + [xl(1);xl(2);xl(1);xl(2)]*plotstruct.basisvectors.vector1(2) + [yl(1);yl(1);yl(2);yl(2)]*plotstruct.basisvectors.vector2(2);
                L = plotstruct.basisvectors.origin(3) + [xl(1);xl(2);xl(1);xl(2)]*plotstruct.basisvectors.vector1(3) + [yl(1);yl(1);yl(2);yl(2)]*plotstruct.basisvectors.vector2(3);
                if isfield(plotstruct.basisvectors,'vector3') && strcmpi(plotstruct.basisvectors.type3,'hkl')
                    H = [H + zl(1)*plotstruct.basisvectors.vector3(1); H + zl(2)*plotstruct.basisvectors.vector3(1)];
                    K = [K + zl(1)*plotstruct.basisvectors.vector3(2); K + zl(2)*plotstruct.basisvectors.vector3(2)];
                    L = [L + zl(1)*plotstruct.basisvectors.vector3(3); L + zl(2)*plotstruct.basisvectors.vector3(3)];
                end
            case 'QEPLANE'
                return; % Not possible!!
            case 'ANGLES'
                return; % Not possible!!
        end
        
        % Create field if not yet existing
        if ~isfield(plotstruct,'linedef'), plotstruct.linedef =  {}; end
        
        % Loop over all sclices
        [ndata, ~, ~, ~, vertinds] = slicecounting(plotstruct.datalist);
        for nsl = 1:ndata
            
            % Get coordinates of plot vertices in HKL
            if strcmpi(plotstruct.type,'hklvectors')
                hkl = plotstruct.vertexlist(vertinds{nsl},1)*plotstruct.basisvectors.vector1 + plotstruct.vertexlist(vertinds{nsl},2) * plotstruct.basisvectors.vector2;
                if isfield(plotstruct.basisvectors,'vector3') && strcmpi(plotstruct.basisvectors.type3,'hkl')
                    hkl = hkl + plotstruct.vertexlist(vertinds{nsl},3) * plotstruct.basisvectors.vector3;
                end
                vertH = hkl(:,1)+plotstruct.basisvectors.origin(1); vertK = hkl(:,2)+plotstruct.basisvectors.origin(2); vertL = hkl(:,3)+plotstruct.basisvectors.origin(3);
            else        
                vl= plotstruct.vertexlist(vertinds{nsl},:);
                if any(strcmpi(plotstruct.type,{'QXY','QXQYEN'}))
                    vl(:,3) = plotstruct.datalist{1}.QVERT;
                end        
                [vertH, vertK, vertL] = calcHKL( vl(:,1), vl(:,2), vl(:,3), UB );
            end

            % Loop over all H,K,L, and create line definition structure hkllindef
            hkllinedef.points = [];
            hkllinedef.connection = [];
            hkllinedef.properties = {'edgecolor','r','linestyle','-'};
            hkllinedef.label = 'HKL';
            hkllinedef.assigntoslice = nsl;
            for h = ceil(min(H)+.01):floor(max(H)-.01)
                [linevert, lineconn] = createintersection(plotstruct.vertexlist(vertinds{nsl},:), plotstruct.datalist{nsl}.faces, vertH - h, 2);
                hkllinedef.connection = [hkllinedef.connection; lineconn + size(hkllinedef.points,1)];
                hkllinedef.points = [hkllinedef.points; linevert];
            end
            for k = ceil(min(K)+.01):floor(max(K)-.01)
                [linevert, lineconn] = createintersection(plotstruct.vertexlist(vertinds{nsl},:), plotstruct.datalist{nsl}.faces, vertK - k, 2);
                hkllinedef.connection = [hkllinedef.connection; lineconn + size(hkllinedef.points,1)];
                hkllinedef.points = [hkllinedef.points; linevert];
            end
            for l = ceil(min(L)+.01):floor(max(L)-.01)
                [linevert, lineconn] = createintersection(plotstruct.vertexlist(vertinds{nsl},:), plotstruct.datalist{nsl}.faces, vertL - l, 2);
                hkllinedef.connection = [hkllinedef.connection; lineconn + size(hkllinedef.points,1)];
                hkllinedef.points = [hkllinedef.points; linevert];
            end
            plotstruct.linedef = [plotstruct.linedef, hkllinedef]; % append this linedef do plotstruct's linedef list
        end
            
        guidata(gcf,plotstruct);
        doplot(plotstruct);
    end

%% Callback: Reset line and plane definitions
    function callback_deleteintersect(src,evnt)
        plotstruct = guidata(gcf);
        if isfield(plotstruct,'linedef'), plotstruct = rmfield(plotstruct,'linedef'); end
        if isfield(plotstruct,'planedef'), plotstruct = rmfield(plotstruct,'planedef'); end
        guidata(gcf,plotstruct);
        doplot(plotstruct);
    end

%% Callback: Export to MFit
    function exportmfit(src,evnt) 
        plotstruct = guidata(gcf);
        x = plotstruct.scandata.x;
        y = plotstruct.scandata.y;
        err = plotstruct.scandata.dy;
        ind = isfinite(y) & isfinite(err);
        xlab = get(get(plotstruct.scanaxeshandle,'xlabel'),'string');
        ylab = get(get(plotstruct.scanaxeshandle,'ylabel'),'string');
        try 
            mfitgo(x(ind),y(ind),err(ind),xlab,ylab,plotstruct.scaninfo);
        catch
            fprintf('Error on exporting to MFit. Check Mfit installation.\n');
        end
    end

%% Callbacks: use nfit
    function gaussfit(src,evnt)
        callnfit('gauss');
    end

    function linearfit(src,evnt)
        callnfit('linear');
    end

    function generalnfit(src,evnt)
        callnfit('*');
    end

    function callnfit(what) 
        plotstruct = guidata(gcf);
        x = plotstruct.scandata.x;
        y = plotstruct.scandata.y;
        err = plotstruct.scandata.dy;
        ind = isfinite(y) & isfinite(err);
        xlab = get(get(plotstruct.scanaxeshandle,'xlabel'),'string');
        ylab = get(get(plotstruct.scanaxeshandle,'ylabel'),'string');
        try 
            figure
            errorbar(x(ind),y(ind),err(ind),'ob');
            xlabel(xlab);
            ylabel(ylab);
            nfit(what);
        catch
            fprintf('Error on exporting to nfit. Check nfit installation.\n');
        end
    end

%% Test
% 
% 
%    function testtest(src,evnt)
%      % Callback for Mouseclick in the 2D-Plot; starts defining a projection line
%      plotstruct = guidata(gcf);
%      [nface,nslice] = findface3d (gca, plotstruct.coordlist, plotstruct.vertexlist, plotstruct.datalist);
%      fprintf('nfac=%d,  nsl=%d\n',nface,nslice);
%    end
% 
% 
%         set(get(plotstruct.axeshandle,'Children'),'ButtonDownFcn',@testtest);

%% End
    if nargout==0, clear fh; end
end %fcplot