function erg = extractsubset(data,field,value,chname)

% function erg = extractsubset(data,field,value)
% Extract a subset from data, defined by data.(field)==val
% For instance: a pal-set from polarized "data"
% or a certain ROI or CHAN from multichannel data
% Guess, at the same time, a plot style by analyzing the pal-Definition
%
% P. Steffens, 11/2016

erg = data;
ind = (data.(field) == value);
erg.coordlist = data.coordlist(ind,:);
erg.valuelist = data.valuelist(ind,:);
if data.raw, erg.monitorlist = data.monitorlist(ind,:); end
if isfield(data,'taglist')
    ii = find(ind); for i = 1:numel(ii), erg.taglist{i} = data.taglist{ii(i)}; end
end
if isfield(data,'pallist'), erg.pallist = data.pallist(ind); end
if isfield(data,'channellist'), erg.channellist = data.channellist(ind); end


% Assign a standard Legend text
if strcmpi(field,'pallist')
    erg.legend = ['pal', num2str(value,'%d'),': '];
elseif strcmpi(field,'channellist')
    if nargin>3, erg.legend = [chname ' ', num2str(value,'%d') ' '];
    else erg.legend = ['Channel ', num2str(value,'%d') ' ']; end
else
    erg.legend = [num2str(value,'%d'),': '];
end

% If polarized data, analyze pal-definition to guess a point style...
if strcmpi(field,'pallist')
    pdfields = fieldnames(data.paldeflist{value});
    for fn=1:length(pdfields)
        val = data.paldeflist{value}.(pdfields{fn});
        if isnumeric(val), val = num2str(val,'%g '); end
        erg.legend = [erg.legend, pdfields{fn}, ' ', val, ' '];
    end

    % Guess a color and point style...
    % (by looking in which direction points H/P1 and looking at flipper states)
    erg.plotstyle.color = 'k'; 
    erg.plotstyle.Marker = 'o';
    erg.plotstyle.MarkerFaceColor = 'none';
    cstring = 'rgb';
    if isfield(data.paldeflist{value},'HX')
        [~,mi] = max(abs(data.paldeflist{value}.HX));
        erg.plotstyle.color = cstring(mi);
    elseif isfield(data.paldeflist{value},'P1')
        [~,mi] = max(abs(data.paldeflist{value}.P1));
        erg.plotstyle.color = cstring(mi);
    elseif isfield(data.paldeflist{value},'fcu')
        if strcmpi(data.paldeflist{value}.fcu,'up'), erg.plotstyle.color = cstring(1); else erg.plotstyle.color = cstring(3); end
    else
        cstring='rkcmgb'; erg.plotstyle.color = cstring(mod(value,6)+1);
    end
    if isfield(data.paldeflist{value},'FLIPPER')
        if strcmpi(data.paldeflist{value}.FLIPPER,'on') || data.paldeflist{value}.FLIPPER==1, erg.plotstyle.MarkerFaceColor = erg.plotstyle.color; end
    elseif ~isempty(strfind(erg.legend,'F2 on')) || ~isempty(strfind(erg.legend,'FLIPPER on')) || ~isempty(strfind(erg.legend,'I10 -'))
        erg.plotstyle.MarkerFaceColor = erg.plotstyle.color; 
    end
    if ~isempty(strfind(erg.legend,'F1 on')) || ~isempty(strfind(erg.legend,'I9 -')), erg.plotstyle.Marker = 's'; end
end
    
erg.type = data.type;
erg.coordtype = data.coordtype;
erg.raw = data.raw;
if isfield(data,'polarized'),    erg.polarized = data.polarized; end
if isfield(data,'multichannel'), erg.multichannel = data.multichannel; end



    