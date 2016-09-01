function [vertexlist, faces, triang, vol] = makevoronoi(points,maxsize,ratio,switches)

% [vertexlist, faces, triang, vol] = makevoronoi(points,maxsize,ratio,switches)
%
% Calculates the voronoi cells for N-dim data points. The maximum size of the
% cells is maxsize in each dimension
% ratio scales the coordinate axes (metrics)
% switches: may contain 'calcvol' and/or 'calctri' to compute volume / delaunay triangulation
%
% vertexlist:    rows are the coordinates of all vertices
% faces:         rows are lists of the vertices (in vertexlist) defining each cell
% triang:        Delaunay triangulation
% vol:           volume / surface area of the cells

% replaces calcpatch.m as generalization to higher dimensions
%
% P. Steffens, 07/2010
%

ndims = size(points,2);     % Number of dimensions
npoints = size(points,1);   % Number of (real) points

calctriang  = (nargin>3) && any(strcmpi(switches,'calctri'));    % Calculate Delaunay triangulation?
calcvol     = (nargin>3) && any(strcmpi(switches,'calcvol'));    % Calculate Volume of cells?


% for each dimension, add two artificial points to make all cells finite
points2 = points;
meanpoint = mean(points,1);
for d = 1:ndims
    s = max(points(:,d)) - min(points(:,d));
    points2 = [points2; meanpoint(1:(d-1)), max(points(:,d)) + s + 3*maxsize(d), meanpoint((d+1):end)];
    points2 = [points2; meanpoint(1:(d-1)), min(points(:,d)) - s - 3*maxsize(d), meanpoint((d+1):end)];
end

% Scale axes
if nargin<3, ratio = ones(1,ndims); end
for d = 1:ndims
    points2(:,d) = points2(:,d) * ratio(d);
end

%compute Voronoi vertices V and cells C
[V,C] = voronoin(points2,{'Qz'}); % ,{'Qbb','Qz'});

% Scale back
for d = 1:ndims
    V(:,d) = V(:,d) / ratio(d);
end





% count vertices to preallocate memory for vertexlist and faces
% (however, not really predictable due to cutting below)
% Note also that this largely overcounted, because it corresponds to a
% non-optimized list (multiple occurences), before "optimizepatch"
nvert=0;
for i=1:length(C)
    nvert = nvert + numel(C{i});
end
% nvert = size(V,1) - 1;  % ("-1" because first point is "inf")
vertexlist=zeros(nvert,ndims); 
faces = zeros(npoints,4);
vol = zeros(npoints,1);


% Obtain Triangulation, if desired
if calctriang
    if ndims == 2
        triang = delaunayfromvoronoi({C{1:npoints}}, points, 2, nvert);
    else
        fprintf('Cannot do triangulation so far for dimension higher than 2!\n'); return;
    end
end


nvert = 1;
for i=1:npoints    
    cut = false;
    newcell = V(C{i},:);
    % Find out if cell needs to be cut - do this for each dimension
    for d=1:ndims
        if max(newcell(:,d)) > points(i,d) + maxsize(d)/2
            if ndims==2, K = [(1:size(newcell,1))', [(2:size(newcell,1))';1]]; else K=[]; end
            [cutpoints, newcell] = polyhypercut( newcell,  [zeros(1,d-1), 1, zeros(1,ndims-d)],   points(i,d) + maxsize(d)/2,    K );
            cut= true;
        end
        if min(newcell(:,d)) < points(i,d) - maxsize(d)/2
            if ndims==2, K = [(1:size(newcell,1))', [(2:size(newcell,1))';1]]; else K=[]; end
            [cutpoints, newcell] = polyhypercut( newcell, -[zeros(1,d-1), 1, zeros(1,ndims-d)], -(points(i,d) - maxsize(d)/2),   K );
            cut = true;
        end
    end
    
    if cut && (ndims ~= 2) %optimize and get Volume
        %(treat ndims=2 below)
        [K,vol(i)] = convhulln(newcell);
        ind = unique(K);
        newcell = newcell(ind,:);
    end
    
    
    % Append vertices of this cell to vertexlist
    vertexlist( nvert : (nvert+size(newcell,1)-1), 1:ndims)  = newcell;
    nvert = nvert + size(newcell,1);
    
    % Append indices of these vertices to cell list (undefined order, not safe for 2D!)
    faces(i,1:size(newcell,1))= (nvert-size(newcell,1)) : (nvert-1);
    
end

faces(faces==0)=NaN;


[vertexlist, faces] = optimizepatch(vertexlist, faces);

if ndims==2 && nargout>2
    vol = abs(polygonarea( vertexlist(:,1), vertexlist(:,2), faces));
end
