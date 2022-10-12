function [zi,errzi] = linearinterpolation(datalist,posi,maxsize)

% Triangle-based linear interpolation in N dimensions of the data in
% datalist at the points posi
%
% Uses Delaunay triangulation
% Slightly modified MATLAB code from "griddatan":
% Uses additionally the errorbars of the given points x,y as wheigt upon averaging 
% maxsize (optional) determines max. size of delaunay triangles in the respective dimension (standard from options file if not given)
%
% P. Steffens, 08/2008


npoints = size(posi,1);
ndims   = size(posi,2);
zi    = NaN(npoints,1); 
errzi = NaN(npoints,1);

z  = datalist.valuelist(:,1); 
dz = datalist.valuelist(:,2); 


if isfield(datalist,'delaunaytri')
    tri = datalist.delaunaytri;
else
    % Triangularize the data
%     tri = delaunayn(datalist.coordlist);      %**
    tri = delaunayfromvoronoi(datalist.faces, datalist.coordlist, 2);
end
if isempty(tri),  warning('linearinterpolation:CannotTriangulate','Data cannot be triangulated.');  return; end
    
% Find those triangles that are too large (unphysically large)
stdcell = getoption(['stdcell.' upper(datalist.coordtype)]);
if nargin < 3
    maxsize = stdcell;
end
toolarge = false(size(tri,1),1);
for d = 1:ndims
    toolarge = toolarge | ( max(reshape(datalist.coordlist(tri,d),size(tri)),[],2) - min(reshape(datalist.coordlist(tri,d),size(tri)),[],2) > maxsize(d) );
end

% Find the nearest triangle (t) and barycentric coordinates (p)
[t,p] = tsearchn(datalist.coordlist,tri,posi);

badts  = ismember(t, find(toolarge))  ; % logical index into t for those triangles that are too large 

% Weighted mean:
for i = 1:npoints
    if ~badts(i) && ~isnan(t(i))
        if datalist.raw == 1
            mo = datalist.monitorlist(:,2);
            [zi(i), errzi(i)] = weightedmean ( z(tri(t(i),:)),  dz(tri(t(i),:)),  p(i,:) .* mo(tri(t(i),:))' );
        else
            [zi(i), errzi(i)] = weightedmean ( z(tri(t(i),:)),  dz(tri(t(i),:)),  p(i,:) ./ dz(tri(t(i),:))'.^2 );
        end
    end
end

if sum(badts) > .1 * npoints     % A significant portion (>10%) of points contain NaN because triangle has been rejected
    fprintf('During interpolation, some points in the result are NaN because data points are too distant.\n');
    if nargin < 3
        fprintf(['If this seems unexpected to you, you may try to increase the value of stdcell.' (upper(datalist.coordtype)) ' in the options file.\n']);
    end
end
