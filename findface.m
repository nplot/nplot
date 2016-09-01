function nface = findface (x, y, coordlist, faces, vertices, delaunaytri)

% Test if point(x,y) lies in any of the faces
% Return NaN or number of face
%
% P. Steffens, 05/2008

nface = NaN;

% If delaunaytri is given, look only at faces around the corner of the
% enclosing triangle
if nargin >5
%     coordlist = coordlist( abs(coordlist(:,1) - x) <= range(1), :);
%     coordlist = coordlist( abs(coordlist(:,2) - y) <= range(2), :);
    tri = tsearch(coordlist(:,1),coordlist(:,2),delaunaytri,x,y);
    if isnan(tri), return; end
    whichfaces = delaunaytri(tri,:);
else
    whichfaces = 1:size(coordlist,1);
end


% Determine distances to points in coordlist
dist = zeros(numel(whichfaces),1);
for nv = 1:numel(whichfaces)
    dist(nv) = sum( (coordlist(whichfaces(nv),:) - [x,y]).^2 );
end

[b,sortedorder] = sort(dist);

% Test faces in order of incresing distance (cannot take only nearest
% because patch is usuall not a Voronoi diagram in qx.qy-coordinates)
for nvmin = whichfaces(sortedorder')
    fac = faces(nvmin, isfinite(faces(nvmin,:)));
    fac = [fac, fac(1)];
    if inpolygon(x,y,vertices(fac,1),vertices(fac,2))
        nface = nvmin;
        return;
    end
end
