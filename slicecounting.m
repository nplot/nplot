function [ndata, npoints, nverts, pointinds, vertinds] = slicecounting(datalist)

% Helper function for indexing points and vertices in plotstruct. Returns:
% ndata:     number of sclices ( length(datalist))
% npoints:   number of points in each slice
% nverts:    number of vertices in each slice
% pointinds: indices of each slice's points in single combined list
% vertinds:  indices of each slice's vertices in single combined list

% P. Steffens, 08/2014

ndata = length(datalist); 
for sl = 1:ndata
    npoints(sl) = size(datalist{sl}.faces,1);  %#ok<*AGROW>
    nverts(sl) = size(datalist{sl}.vertexlist,1);
    pointinds{sl} = (1:npoints(sl)) + sum(npoints(1:sl-1));
    vertinds{sl} = (1:nverts(sl)) + sum(nverts(1:sl-1));
end 