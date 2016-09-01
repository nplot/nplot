function plot1d(list,axdim,axhandle,varargin)

% function plot1d(list,axdim,axhandle)
% list:     Data in liststruct format, may be cell array
% axdim:    Which coordinate dimenstion to use as x-axis
% axhandle: given in case previous axes are to be used 
%           (can be array to use multiple when list is cell array)

% P. Steffens, 06/2012

if ~iscell(list), m = list; clear list; list{1}=m;  clear m; end % ensure cell array

plottype={'ob','or','ok','og','oc','fb','fr','fk','fg','fc'};
multipleaxes = false;
yoffset = readinput('offset',varargin); if isempty(yoffset), yoffset = 0; end

if nargin<3 || isempty(axhandle)
    figure
else 
    if numel(axhandle)>1
        if numel(axhandle) == length(list), multipleaxes = true;
        else fprintf('Error: Number of axes handles not consistent with number of plots.\n'); return; end
        if ~all(ishandle(axhandle)), fprintf('Error: invalid axes.\n'); return; end
    end
    set(gcf,'currentaxes',axhandle(1));    
end



hold on
box on



for num = 1:length(list)

    if multipleaxes
        set(gcf,'currentaxes',axhandle(num));
        hold on
        box on
    end
    
    figdata = guidata(gca);
    if isempty(figdata) || ~all(isfield(figdata,{'clist','mlist','tlist','datalegendtext','datahandles','fitlinehandles'}))
        figdata.clist=[];
        figdata.mlist=[];
        figdata.tlist={};
        figdata.datalegendtext = {};
        figdata.datahandles = [];
        figdata.fitlinehandles = [];
        skiplegend=0;
        numdatsets = 0;
    else    
        numdatsets = numel(get(gca,'children'));
    end

    filled = false;

    if isempty(list{num}) || isempty(list{num}.coordlist), continue; end
    if isfield(list{num},'plotstyle') && ischar(list{num}.plotstyle)
         linespec = list{num}.plotstyle; 
    else linespec = plottype{mod(num-1+numdatsets,length(plottype))+1}; 
    end
    if strfind(linespec,'f'), linespec(linespec=='f') = 'o'; filled = true; else filled=false; end
    h = errorbar(list{num}.coordlist(:,axdim), list{num}.valuelist(:,1)+yoffset, list{num}.valuelist(:,2), linespec , ...
                'MarkerFaceColor', 'auto', 'MarkerSize',8, 'Linewidth',1.3);
    if filled, set(h,'MarkerfaceColor',get(h,'Color')); else set(h,'MarkerfaceColor','none'); end
    set(h,'tag','data');
    set(h,'ButtonDownFcn',@mouseclick1); 
    figdata.datahandles = [figdata.datahandles, h];
    if isfield(list{num},'plotstyle') && isstruct(list{num}.plotstyle) % set custom properties
        pfields = fieldnames(list{num}.plotstyle);
        for fn=1:length(pfields)
            set(h, pfields{fn}, list{num}.plotstyle.(pfields{fn}));
        end
    end
    figdata.clist = [figdata.clist; list{num}.coordlist(:,axdim), list{num}.valuelist(:,1)];
    if isfield(list{num},'monitorlist'),
        figdata.mlist = [figdata.mlist; list{num}.monitorlist(:,1)];
    else
        figdata.mlist = [figdata.mlist; zeros(size(list{num}.coordlist,1),1)];
    end
    if isfield(list{num},'taglist'),
        figdata.tlist = [figdata.tlist, list{num}.taglist{:}];
    else
        emptytag = cell(size(list{num}.coordlist,1),1);
        figdata.tlist = {figdata.tlist{:}, emptytag{:}};
    end
    if isfield(list{num},'legend'),
        figdata.datalegendtext = {figdata.datalegendtext{:}, list{num}.legend};
    else
        figdata.datalegendtext = {figdata.datalegendtext{:}, []};
    end

    
    if ~any(strcmpi(varargin,'nolegend'))   % Show legend
        leg = {}; lh = [];
        for li=1:numel(figdata.datahandles)     
            if ~ isempty(figdata.datalegendtext{li})
                lh = [lh, figdata.datahandles(li)];
                leg = {leg{:}, figdata.datalegendtext{li}};
            end
        end        
        legend(lh, leg);
    end

    set(gca,'ButtonDownFcn',@mouseclick1)

    guidata(gca, figdata);

end

%% Callback routine for Mouseclick in Axes

   function mouseclick1(src,evnt)
     % Callback for Mouseclick 
     if strcmp(get(gcf,'SelectionType'),'normal')
        figdata = guidata(gca);
        cp = get(gca,'CurrentPoint');
        xmouse = cp(1,1); ymouse = cp(1,2);
        xl = xlim; yl = ylim; xr = xl(2)-xl(1); yr = yl(2)-yl(1);
        [m,npoint] = min(((figdata.clist(:,1)-xmouse)/xr).^2 + ((figdata.clist(:,2)-ymouse)/yr).^2);
        fprintf(['x=' num2str(figdata.clist(npoint,1)) ', y=' num2str(figdata.clist(npoint,2))]);
        if figdata.mlist(npoint)>0
            fprintf([', Mon=' num2str(figdata.mlist(npoint))]);
        end
        if ~isempty(figdata.tlist(npoint))
            fprintf([', Scans # ' figdata.tlist{npoint}]);
        end
        fprintf('\n');
     end
   end


end

        