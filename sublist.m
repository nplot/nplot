function erg = sublist(list, points)

% erg = sublist(list, points)
% Return a subset (defined by points) of the list structure

% P. Steffens, 09/2008 - 09/2014


list.coordlist = list.coordlist(points,:);
list.valuelist = list.valuelist(points,:);

if isfield(list,'monitorlist')
    list.monitorlist = list.monitorlist(points,:);
end

% Remove unused vertices
if all(hasfield(list,{'vertexlist','faces'}))
    list.faces     = list.faces(points,:);
    [list.vertexlist, list.faces] = optimizepatch( list.vertexlist, list.faces);
end

if isfield(list,'delaunaytri')
    list = rmfield(list,'delaunaytri');
end

if isfield(list,'taglist')
    list.taglist = {list.taglist{points}};
end

if isfield(list,'pallist')
    list.pallist = list.pallist(points,:);
end

if isfield(list,'channellist')
    list.channellist = list.channellist(points,:);
end

erg = list;