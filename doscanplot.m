function doscanplot(plotstruct)

% Create and plot 1D scan

% P. Steffens, 03/2008
% Commented parts are now in makescandata.m


% %%
% if ~isfield(plotstruct.scandef,'xdat'), return; end
% startpoint = [plotstruct.scandef.xdat(1), plotstruct.scandef.ydat(1)];
% endpoint   = [plotstruct.scandef.xdat(2), plotstruct.scandef.ydat(2)];
% np         =  plotstruct.scandef.np;
% 
% %% construct the scan points
% scanpath = repmat(startpoint,np,1) + (0.5:(np-0.5))'/np  * (endpoint-startpoint);
% 
% %% Obtain values and errors by interpolation, integration or pojection
% 
% 
% if strcmpi(plotstruct.scandef.type,'Interpolation')
%     datalist = plotstruct.datalist;
%     datalist.coordlist = plotstruct.coordlist;
%     datalist.coordtype = plotstruct.type; % ** Translations... unfortunately... conventions should be changed...
%         if strcmpi(plotstruct.type,'qeplane'), datalist.coordtype = 'A4ENERGY'; end
%         if strcmpi(plotstruct.type,'qxqyen'),  datalist.coordtype = 'ANGLESENERGY'; end
%     [plotstruct.scandata.y, plotstruct.scandata.dy]  = linearinterpolation (datalist, [scanpath(:,1), scanpath(:,2)]);
% end
% 
% if strcmpi(plotstruct.scandef.type,'Integration')
%     integ = integratepatch (plotstruct.datalist.faces, plotstruct.vertexlist, plotstruct.datalist.valuelist(:,1), plotstruct.datalist.valuelist(:,2), ...
%                             startpoint - plotstruct.scandef.framevec, [endpoint-startpoint; 2 * plotstruct.scandef.framevec], [np,1] );
%     sd  = integ{1};
%     plotstruct.scandata.y = sd(:,1);  plotstruct.scandata.dy = sd(:,2);                   
% end
% 
% if strcmpi(plotstruct.scandef.type,'Projection')
%     % Find all points within the parallelogramm
%     % Now cut along the two edges parallel to framevec:
%     normal = [-plotstruct.scandef.framevec(2), plotstruct.scandef.framevec(1)]; normal = normal' ./ sqrt(normal*normal');
%     C1 = startpoint * normal;
%     C2 = endpoint   * normal;
%     % Indizes of those points between the two edges
%     between1 = (plotstruct.coordlist(:,1) * normal(1) + plotstruct.coordlist(:,2) * normal(2) < C1) == ...
%                (plotstruct.coordlist(:,1) * normal(1) + plotstruct.coordlist(:,2) * normal(2) > C2);
%     % Cut along the two edges parallel to the scan line:
%     normal = endpoint-startpoint; normal = [-normal(2), normal(1)]; normal = normal' ./ sqrt(normal*normal');
%     C1 = (startpoint + plotstruct.scandef.framevec) * normal;
%     C2 = (startpoint - plotstruct.scandef.framevec) * normal;
%     % Indizes of those points between the two edges
%     between2 = (plotstruct.coordlist(:,1) * normal(1) + plotstruct.coordlist(:,2) * normal(2) < C1) == ...
%                (plotstruct.coordlist(:,1) * normal(1) + plotstruct.coordlist(:,2) * normal(2) > C2);
% 
%     % These points make the scan
%     plotstruct.scandata.y  = plotstruct.datalist.valuelist(between1 & between2, 1);
%     plotstruct.scandata.dy = plotstruct.datalist.valuelist(between1 & between2, 2);     
%     
%     % Find the projection (along framevec) of these points on the scanline
%     C = startpoint * normal;
%     ind = find(between1 & between2);
%     scanpath = zeros(0,2);
%     for i=1:numel(ind)
%         scanpath(i,:) = plotstruct.coordlist(ind(i),:) + (C - plotstruct.coordlist(ind(i),:)*normal)/(plotstruct.scandef.framevec*normal) * plotstruct.scandef.framevec;
%     end    
% end

%% Set current axes
axes(plotstruct.scanaxeshandle);
cla(plotstruct.scanaxeshandle);
hold on
box on

% if isempty(scanpath), return; end

% %% Determine x-coordinate
% 
% % first, calculate HKL values
% UB = UBmatrix( plotstruct.datalist.sampleinfo.lattice, plotstruct.datalist.sampleinfo.ax, plotstruct.datalist.sampleinfo.bx);
% [H, K, L] = calcHKL( scanpath(:,1), scanpath(:,2), -plotstruct.datalist.QVERT, UB );  % ** '-'
% 
% % Now, x-coordinate depending on setting of scandef.xaxiscoord:
% switch plotstruct.scandef.xaxiscoord
%     case 'AUTO'
%         [m,i] = max( max(scanpath,[],1) - min(scanpath,[],1)); %i is No. of column with largest variation; take this as x-axis
%         xlabels=[];
%         if strcmpi(plotstruct.type, 'qxy'),  
%             xlabels= {['Q_x (' char(197) '^{-1})'],['Q_y (' char(197) '^{-1})']};
%             plotstruct.scandata.x = scanpath(:,i);
%             if getoption('plotopt.preferHKL') % Do the plot in HKL coords instead Angstroms
%                 HKL = [H, K, L];
%                 [m,i] = max( max(HKL,[],1) - min(HKL,[],1)); %i is No. of column with largest variation; take this as x-axis
%                 xlabels= {'H (r.l.u)', 'K (r.l.u)', 'L (r.l.u)'};
%                 plotstruct.scandata.x = HKL(:,i);
%             end
%         elseif strcmpi(plotstruct.type, 'qeplane')
%             xlabels={['|Q| (' char(197) '^{-1})'], 'Energy (meV)'}; 
%             plotstruct.scandata.x = scanpath(:,i);
%         elseif strcmpi(plotstruct.type, 'angles')
%             xlabels={'Scattering angle (in plane)', 'Sample rotation angle'}; 
%             plotstruct.scandata.x = scanpath(:,i);
%         else
%             plotstruct.scandata.x = scanpath(:,i);
%         end
%         if ~isempty(xlabels), xlabel(xlabels{i}); end
%         
%     case 'QX'
%         xlabel(['Q_x (' char(197) '^{-1})']);
%         plotstruct.scandata.x = scanpath(:,1);
%     case 'QY'
%         xlabel(['Q_y (' char(197) '^{-1})']);
%         plotstruct.scandata.x = scanpath(:,2);
%     case 'QH'
%         xlabel('H (r.l.u)');
%         plotstruct.scandata.x = H;
%     case 'QK'
%         xlabel('K (r.l.u)');
%         plotstruct.scandata.x = K;
%     case 'QL'
%         xlabel('L (r.l.u)');
%         plotstruct.scandata.x = L;
%     case 'QM'
%         xlabel(['|Q| (' char(197) '^{-1})']);
%         plotstruct.scandata.x = scanpath(:,1);
%     case 'EN'
%         xlabel('Energy (meV)');
%         plotstruct.scandata.x = scanpath(:,2);
% end

[plotstruct.scandata, xlabeltext] = makescandata(plotstruct);


%% Plot
errorbar(plotstruct.scandata.x, plotstruct.scandata.y, plotstruct.scandata.dy,'ob');
if isfield(plotstruct,'properties') && isfield(plotstruct.properties,'normalization')
    ylabel(['Counts per ' plotstruct.properties.normalization]);
end
xlabel(xlabeltext);

%%
guidata(gcf, plotstruct);