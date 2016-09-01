function  clickbgr (xx,yy,dy)



plotdata.xpoint = [xx(1), xx(end)];
plotdata.ypoint = [yy(1), yy(end)];
plotdata.xx = xx;
plotdata.yy = yy;
plotdata.dy = dy;


figure
plotdata.axes1handle = subplot(2,1,1);
hold on

errorbar(plotdata.xx, plotdata.yy, plotdata.dy, 'ob');

get(plotdata.axes1handle, 'ylim');

plotdata.plinex = linspace(min(plotdata.xpoint), max(plotdata.xpoint), 1000);
plotdata.pliney = spline(plotdata.xpoint, plotdata.ypoint, plotdata.plinex);        
plotdata.plineh = plot(plotdata.plinex, plotdata.pliney, '-r');
plotdata.ppointsh = plot(plotdata.xpoint, plotdata.ypoint, 'or');

set(plotdata.axes1handle, 'YLimMode', 'manual');        

set(plotdata.axes1handle,'ButtonDownFcn',@addpoint);
set(get(plotdata.axes1handle,'Children'),'ButtonDownFcn',@addpoint);
set(plotdata.ppointsh,'ButtonDownFcn',@movepoint);

plotdata.axes2handle = subplot(2,1,2);

showdifference;

guidata(gcf,plotdata);

%createline;

%%

%     function createline
%         plotdata = guidata(gcf);
%         
%         axes(plotdata.axes1handle);
%         hold off
%         plotdata.ylim = get(plotdata.axes1handle, 'ylim');
%         errorbar(plotdata.xx, plotdata.yy, plotdata.dy, 'ob');
%         
%         plotdata.plineh = plot(plotdata.plinex, plotdata.pliney, '-r');
%         plotdata.ppointsh = plot(plotdata.xpoint, plotdata.ypoint, 'or');
%         set(plotdata.axes1handle, 'ylim', plotdata.ylim);
%         
%         set(plotdata.ppointsh,'ButtonDownFcn',@movepoint);
%         
%         bgr = spline(plotdata.xpoint, plotdata.ypoint, xx);
%         
%         showdifference;
%         guidata(gcf,plotdata);
%     end

    function showdifference
        axes(plotdata.axes2handle);
        plotdata.diff = plotdata.yy - spline(plotdata.xpoint, plotdata.ypoint, plotdata.xx);
        errorbar(plotdata.xx, plotdata.diff, plotdata.dy, 'ob');
    end
        



    function addpoint(src,evnt)
    % Callback for Mouseclick in the 2D-Plot; adds a point
         if strcmpi(get(gcf,'SelectionType'),'normal') %Add   
            plotdata = guidata(gcf);
            cp = get(plotdata.axes1handle,'CurrentPoint');
          	xp = cp(1,1); yp = cp(1,2);        
            plotdata.xpoint = [plotdata.xpoint, xp];
            plotdata.ypoint = [plotdata.ypoint, yp];        
            % Calculate new values of BG for the line
            plotdata.pliney = spline(plotdata.xpoint, plotdata.ypoint, plotdata.plinex);
            % Update the plot
            set(plotdata.ppointsh, 'XData', plotdata.xpoint, 'YData', plotdata.ypoint);
            set(plotdata.plineh, 'YData', plotdata.pliney); drawnow expose
            showdifference;
            guidata(gcf,plotdata);
         end
    end
    

    function movepoint(src,evnt)
        % Callback for moving or deleting an existing point
        cp = get(plotdata.axes1handle,'CurrentPoint'); 
        xp = cp(1,1); yp = cp(1,2);
        if strcmp(get(gcf,'SelectionType'),'normal')  % Move        
            % Identify the point to be moved
            [m,pointnr] = min( (plotdata.xpoint - xp).^2 + (plotdata.ypoint - yp).^2 );
            set(gcf,'WindowButtonMotionFcn',@mousemove1)
            set(gcf,'WindowButtonUpFcn',@mouseup1)
         elseif strcmpi(get(gcf,'SelectionType'),'alt') % Delete
            % Identify the point to be moved
            [m,pointnr] = min( (plotdata.xpoint - xp).^2 + (plotdata.ypoint - yp).^2 );
            plotdata.xpoint = plotdata.xpoint([1:pointnr-1, pointnr+1:end]);
            plotdata.ypoint = plotdata.ypoint([1:pointnr-1, pointnr+1:end]); 
             % Calculate new values of BG for the line
            plotdata.pliney = spline(plotdata.xpoint, plotdata.ypoint, plotdata.plinex);
            % Update the plot
            set(plotdata.ppointsh, 'XData', plotdata.xpoint, 'YData', plotdata.ypoint);
            set(plotdata.plineh, 'YData', plotdata.pliney); drawnow expose
            showdifference;
            guidata(gcf,plotdata);
        end
 
        function mousemove1(src,evnt)
           cp = get(plotdata.axes1handle,'CurrentPoint');
           plotdata.xpoint(pointnr) = cp(1,1);
           plotdata.ypoint(pointnr) = cp(1,2);
           % Calculate new values of BG for the line
           plotdata.pliney = spline(plotdata.xpoint, plotdata.ypoint, plotdata.plinex);
           % Update the plot
           set(plotdata.ppointsh, 'XData', plotdata.xpoint, 'YData', plotdata.ypoint);
           set(plotdata.plineh, 'YData', plotdata.pliney); drawnow expose
           % Update difference plot           
           showdifference;
        end
   
        function mouseup1(src,evnt)
           % Reset behaviour of figure
           set(src,'WindowButtonMotionFcn','')
           set(src,'WindowButtonUpFcn','')
           % Calculate new values of BG at the data x values
           plotdata.ybgr = spline(plotdata.xpoint, plotdata.ypoint, plotdata.xx);
           guidata(gcf, plotdata);
        end
    end

end
