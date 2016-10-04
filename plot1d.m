function linespecs =  plot1d(list,varargin)

% function linespecs = plot1d(list,axdim,axhandle)
% list:     Data in liststruct format, may be cell array
% varargin may contain:
% 'axdim', dim:         Which coordinate dimension to use as x-axis
% 'axhandle', handle:   given in case previous axes are to be used 
%           (can be array to use multiple when list is cell array)
% 'offset', yoff:       Offset on y-axis
% 'plotstyle', pstyle:  (default) plotstyle string. (Overrides .plotstyle field in list) 
% Returns:
% linespecs: Cell array of plotstyles

% - saves a struct plotdataset in guidata of current figure for each data set:
%   .x, .y, .dy, .axhandle, .plothandle, (.mlist), (.tlist), (.legendtext)
% - updates the data cursor display of the figure window appropriately

% P. Steffens, 06/2016


if ~iscell(list), m = list; clear list; list{1}=m;  clear m; end % ensure cell array

plottype={'ob','or','ok','og','oc','fb','fr','fk','fg','fc'};
multipleaxes = false;
linespecs = {};

% interpret varargin
axdim   = readinput('axdim',varargin); if isempty(axdim), axdim = 1; end
yoffset = readinput('offset',varargin); if isempty(yoffset), yoffset = 0; end
pstyle  = readinput('plotstyle',varargin);
axhandle = readinput('axhandle',varargin);
if isempty(axhandle)
    figure;     axhandle = axes;
else 
    if numel(axhandle)>1
        if numel(axhandle) == length(list), multipleaxes = true;
        else fprintf('Error: Number of axes handles not consistent with number of plots.\n'); return; end
        if ~all(ishandle(axhandle)), fprintf('Error: invalid axes.\n'); return; end
    end 
end
set(gcf,'currentaxes',axhandle(1));     hold on;    box on

% get guidata from figure
figdata = guidata(gca);
if ~isfield(figdata,'plotdataset'), figdata.plotdataset = {}; end
dcount = length(figdata.plotdataset); 

%Loop over datasets
for num = 1:length(list)
    
    if isempty(list{num}) || isempty(list{num}.coordlist), continue; end
    if multipleaxes, set(gcf,'currentaxes',axhandle(num)); hold on; box on; end % set axis, if necessary
   
    % set plotstyle
    if ~isempty(pstyle), list{num}.plotstyle=pstyle; end
    % use custom plotstyle if string. (Can also be struct, see below)
    if isfield(list{num},'plotstyle') && ischar(list{num}.plotstyle), linespec = list{num}.plotstyle; 
    else linespec = plottype{mod(num-1+dcount,length(plottype))+1}; % default
    end 
    if strfind(linespec,'f'), linespec(linespec=='f') = 'o'; filled = true; else filled=false; end
    
    linespecs{num}=linespec; %#ok<AGROW>
    
    % Plot
    h = errorbar(list{num}.coordlist(:,axdim), list{num}.valuelist(:,1)+yoffset, list{num}.valuelist(:,2), linespec , ...
                'MarkerFaceColor', 'auto', 'MarkerSize',8, 'Linewidth',1.3);
    if filled, set(h,'MarkerfaceColor',get(h,'Color')); else set(h,'MarkerfaceColor','none'); end
    % if plotstyle is struct, set custom properties now
    if isfield(list{num},'plotstyle') && isstruct(list{num}.plotstyle) 
        pfields = fieldnames(list{num}.plotstyle);
        for fn=1:length(pfields)
            set(h, pfields{fn}, list{num}.plotstyle.(pfields{fn}));
        end
    end
    
    set(h,'tag','data'); %(cursorfcn (below) looks for this)
%     if any(strcmpi(varargin,'printcallback')), set([h,gca],'ButtonDownFcn',@mouseclick1); end
    set(datacursormode(gcf),'UpdateFcn',@cursorfcn);
   
    % set up structure for this data set
    thisplotdata.axhandle = axhandle; % "parent" axes
    thisplotdata.plothandle = h; % handle of errorbar object
    thisplotdata.x = list{num}.coordlist(:,axdim);
    thisplotdata.y = list{num}.valuelist(:,1);
    if size(list{num}.valuelist,2)<=2,      thisplotdata.dy = list{num}.valuelist(:,2); end
    if isfield(list{num},'monitorlist'),    thisplotdata.mlist = list{num}.monitorlist(:,1); end
    if isfield(list{num},'taglist'),        thisplotdata.tlist = list{num}.taglist; end
    if isfield(list{num},'legend'),         thisplotdata.legendtext = list{num}.legend; end

    figdata.plotdataset = [ figdata.plotdataset, thisplotdata];
    
    % Show legend
    if ~any(strcmpi(varargin,'nolegend')) && isfield(list{num},'legend')    
        if isempty(legend(gca))
            legend(h, list{num}.legend); 
        else % get current legend and add new entry
            if verLessThan('matlab','8.4') %<R2014b
                [~, ~, legobjs,text_strings]= legend(gca); 
            else % new syntax of Matlab legend command
                ax = gca;
                legobjs = findobj(ax.Children,'flat','-regexp','Displayname','\S+');
                text_strings = {};
                for j=1:numel(legobjs), text_strings{j} = legobjs(j).DisplayName; end %#ok<AGROW>
                legobjs = legobjs(end:-1:1); text_strings = text_strings(end:-1:1);
            end
            legend([legobjs', h], [ text_strings, {list{num}.legend}, ]);
        end
    end
    
end

guidata(gcf, figdata);

%% Update function for data cursor
    function txt = cursorfcn(~,event_obj)
        % Customizes text of data tips
        pos = get(event_obj,'Position');
        txt = {['X: ',num2str(pos(1))], ['Y: ',num2str(pos(2))]}; % default text
        if get(datacursormode(gcf),'snaptodatavertex') % try to get additional info
            tar = get(event_obj,'Target'); % handle of graphics object of data tip
            if ~strcmpi(get(tar,'tag'),'data'), return; end
            figdata = guidata(gcf);
            if ~isfield(figdata,'plotdataset'), return; end
            ind = 1; l = length(figdata.plotdataset);
            while ind<=l && figdata.plotdataset{ind}.plothandle ~= tar, ind=ind+1; end
            if ind>l, return; end % object not found
            % identify index of data point (by min. distance)
            [~,npoint] = min((figdata.plotdataset{ind}.x-pos(1)).^2 + (figdata.plotdataset{ind}.y-pos(2)).^2);
            if isfield(figdata.plotdataset{ind},'dy'),    txt = [txt, {['dY: ',num2str(figdata.plotdataset{ind}.dy(npoint))]}]; end  
            if isfield(figdata.plotdataset{ind},'mlist'), txt = [txt, {['Sum Mon.: ',num2str(figdata.plotdataset{ind}.mlist(npoint))]}]; end  
            if isfield(figdata.plotdataset{ind},'tlist'), txt = [txt, {['Data: ',num2str(figdata.plotdataset{ind}.tlist{npoint})]}]; end  
        end
    end


% %% Callback routine for Mouseclick in Axes
% 
%    function mouseclick1(src,evnt)
%      % Callback for Mouseclick 
%      if strcmp(get(gcf,'SelectionType'),'normal')
%         figdata = guidata(gca);
%         cp = get(gca,'CurrentPoint');
%         xmouse = cp(1,1); ymouse = cp(1,2);
%         xl = xlim; yl = ylim; xr = xl(2)-xl(1); yr = yl(2)-yl(1);
%         [m,npoint] = min(((figdata.clist(:,1)-xmouse)/xr).^2 + ((figdata.clist(:,2)-ymouse)/yr).^2);
%         fprintf(['x=' num2str(figdata.clist(npoint,1)) ', y=' num2str(figdata.clist(npoint,2))]);
%         if figdata.mlist(npoint)>0
%             fprintf([', Mon=' num2str(figdata.mlist(npoint))]);
%         end
%         if ~isempty(figdata.tlist(npoint))
%             fprintf([', Scans # ' figdata.tlist{npoint}]);
%         end
%         fprintf('\n');
%      end
%    end


end

        